# Helper functions
##################

sourceCpp("helper.cpp")

normalize <- function(x) {x/rowSums(x)}

gen_data_village <- function(params_true, W_0, ps = 0, q = 0) {
  N <- params_true$N
  base_prob <- params_true$base_prob
  mu <- params_true$mu
  deg_beta <- params_true$deg_beta
  sigma2 <- params_true$sigma2

  # Randomize treatment and outcomes
  itt <- rbinom(N, 1, base_prob)
  exposure <- itt*2 + (W_0 %*% itt > 0) + 1

  y_0 <- rnorm(N, mu[exposure] + rowSums(W_0)*deg_beta[exposure], sigma2[exposure])

  gen_data <- list()
  for(i in 1:length(ps)) {
    fakeEdges = rbinom(N*N, 1, q)
    lostEdges = rbinom(N*N, 1, ps[i])

    #Observed corrupted network
    W_0_mod <- matrix(data = pmax(c(W_0) - lostEdges,0) +  ((c(W_0) == 0) & (fakeEdges == 1)), nrow = N, ncol = N)
    diag(W_0_mod) = 0

    gen_data[[i]] = list(y_0 = y_0, W_0_true = W_0, W_0 = W_0_mod, itt = itt, itt2 = itt, subset = rep(TRUE,N))
  }
  return(gen_data)
}





# Horvitz Thompson estimator from Aronow and Samii (2017)
#########################################################

est_effects_ht <- function(data, base_prob) {
  W_0 <- data$W_0

  # Restrict to subjects with at least one incoming edge
  hasFriend <- rowSums(W_0) > 0
  W_0 <- W_0[hasFriend, ]
  itt <- data$itt
  y_0 <- data$y_0[hasFriend]
  itt_subset <- itt[hasFriend]

  # Calculate and estimate probability of exposure
  est_exposure <- itt_subset*2 + (W_0 %*% itt > 0) + 1
  degrees <- rowSums(W_0)
  prob_exposure <- matrix(data = 0, nrow = length(y_0), ncol = 4)
  prob_exposure[,1] = (1-base_prob)^(1+degrees)
  prob_exposure[,2] = (1-base_prob)*(1-(1-base_prob)^degrees)
  prob_exposure[,3] = base_prob*(1-base_prob)^(degrees)
  prob_exposure[,4] = base_prob*(1-(1-base_prob)^(degrees))

  mu_est <- sapply(1:4, function(i) sum(y_0[est_exposure == i]/prob_exposure[est_exposure == i,i])/length(y_0))

  return(mu_est)
}





# Proposed EM estimator (linear)
################################

est_effects_em <- function(data, p_prob, q_prob, alpha, beta,
                           N2 = NA, num_treated = integer(0),
                           min_missing = -10, max_missing = 10, max_iter = 200, epsilon = 0.001) {
  subset = data$subset
  N <- sum(subset)
  y <- data$y_0[subset]
  itt <- data$itt[subset]

  # Allow for indirect treatment to depend on a different treatment assignment than direct treatment
  itt2 <- data$itt2
  trt_friends = (data$W_0 %*% itt2)[subset]
  num_friends = rowSums(data$W_0)[subset]
  if(length(num_treated) == 0) {
    num_treated = rep(sum(data$itt2),N)
  }
  if(is.na(N2)) {
    N2 <- length(data$itt2)
  }
  itt2 <- itt2[subset]



  results <- list()

  num_friends_x = sapply(min_missing:max_missing, function(x) num_friends + x)
  xlen = max_missing - min_missing + 1
  old_log_lik = -Inf

  # Estimate prior distribution of exposure conditions x true degrees (lexographical order)
  raw_prior <- t(est_clprior(min_missing, max_missing, itt, itt2, N2, num_treated, trt_friends, num_friends, p_prob, q_prob, alpha, beta))
  clprior_est <- normalize(raw_prior)

  # Initialize parameters based on prior
  params_est = est_params_eqvar(clprior_est, num_friends_x, y)

  for(iter in 1:max_iter) {
    # Estimate posterior distribution of exposure conditions x true degrees
    raw_post <- est_clpost(clprior_est, num_friends_x, y, params_est)
    log_lik <- sum(log(rowSums(raw_post)))
    clpost_est <- normalize(raw_post)

    # Marginal posterior of true degrees
    post_mat = t(clpost_est[,1:xlen] + clpost_est[,xlen + 1:xlen] + clpost_est[,2*xlen + 1:xlen] + clpost_est[,3*xlen + 1:xlen])

    # Estimate corruption in network p and q
    opt = optim(fn = corr_objective, gr = corr_gradient, par = c(p_prob,q_prob),
                post_mat = post_mat, num_friends = num_friends, N2 = N2, min_missing = min_missing, alpha = alpha, beta = beta,
                method = "L-BFGS-B", lower = c(1e-06,1e-06), upper = c(0.9,0.05), control = list(factr = 1e12))
    p_prob = opt$par[1]
    q_prob = opt$par[2]

    # Estimate parameters
    params_est <- est_params_eqvar(clpost_est, num_friends_x, y)

    # Update prior distribution
    raw_prior <- t(est_clprior(min_missing, max_missing, itt, itt2, N2, num_treated, trt_friends, num_friends, p_prob, q_prob, alpha, beta))
    clprior_est <- normalize(raw_prior)



    if(abs(log_lik - old_log_lik) > epsilon) {
      old_log_lik = log_lik
    } else {
      x_range = dim(num_friends_x)[2]
      x_mean0 = sum(num_friends_x * (clpost_est[,1:x_range] + clpost_est[,x_range + 1:x_range] +
                                      clpost_est[,2*x_range + 1:x_range] + clpost_est[,3*x_range + 1:x_range]))/
        sum((num_friends_x >= 0) * (clpost_est[,1:x_range] + clpost_est[,x_range + 1:x_range] +
                                      clpost_est[,2*x_range + 1:x_range] + clpost_est[,3*x_range + 1:x_range]))
      results$means0 = params_est[2,] + x_mean0*params_est[1,]
      results$params = params_est
      results$p_prob = p_prob
      results$q_prob = q_prob
      results$log_lik = log_lik
      results$iter = iter
      results$eff_obs = c(sum((num_friends_x > 0)*clpost_est[,1:x_range]), sum((num_friends_x > 0)*clpost_est[,x_range+1:x_range]),
                          sum((num_friends_x > 0)*clpost_est[,2*x_range+1:x_range]), sum((num_friends_x > 0)*clpost_est[,3*x_range+1:x_range]))
    }
  }

  return(results)
}










# Proposed EM estimator (logistic)
##################################

est_effects_em_binom <- function(data, p_prob, q_prob, alpha, beta,
                           min_missing = -10, max_missing = 10, max_iter = 200, epsilon = 0.001) {
  N <- data$N
  N2 <- data$N2
  num_correct <- data$num_correct
  num_attempt <- data$num_attempt
  num_treated <- data$num_treated
  itt <- data$itt
  itt2 <- data$itt2
  trt_friends <- data$trt_friends
  num_friends <- data$num_friends

  results <- list()

  num_friends_x = sapply(min_missing:max_missing, function(x) num_friends + x)
  xlen = max_missing - min_missing + 1
  old_log_lik = -Inf

  # Estimate prior distribution of exposure conditions x true degrees (lexographical order)
  raw_prior <- t(est_clprior(min_missing, max_missing, itt, itt2, N2, num_treated, trt_friends, num_friends, p_prob, q_prob, alpha, beta))
  clprior_est <- normalize(raw_prior)

  # Initialize parameters based on prior
  params_est <- matrix(data = 0, nrow = 2, ncol = 4)
  for(cl in 0:3) {
    params_est[,cl + 1] = optim(par = params_est[,cl+1], fn = param_binom_objective, post_mat = clprior_est[,1:xlen+cl*xlen],
                                num_friends_x = num_friends_x, n = num_attempt, a = num_correct,
                                method = "BFGS", control = list(factr = 1e12))$par
  }

  for(iter in 1:max_iter) {
    # Estimate posterior distribution of exposure conditions x true degrees
    raw_post <- est_clpost_binom(clprior_est, num_friends_x, n = num_attempt, a = num_correct, params_est)
    log_lik <- sum(log(rowSums(raw_post)))
    clpost_est <- normalize(raw_post)

    # Marginal posterior of true degrees
    post_mat = t(clpost_est[,1:xlen] + clpost_est[,xlen + 1:xlen] + clpost_est[,2*xlen + 1:xlen] + clpost_est[,3*xlen + 1:xlen])

    # Estimate corruption in network p and q
    opt = optim(fn = corr_objective, gr = corr_gradient, par = c(p_prob,q_prob),
                post_mat = post_mat, num_friends = num_friends, N2 = N2, min_missing = min_missing,
                alpha = alpha, beta = beta, method = "L-BFGS-B", lower = c(1e-06,1e-06), upper = c(0.7,0.05), control = list(factr = 1e12))
    p_prob = opt$par[1]
    q_prob = opt$par[2]

    # Estimate parameters
    for(cl in 0:3) {
      params_est[,cl + 1] = optim(par = params_est[,cl+1], fn = param_binom_objective, post_mat = clpost_est[,1:xlen+cl*xlen],
                                  num_friends_x = num_friends_x, n = num_attempt, a = num_correct,
                                  method = "BFGS", control = list(factr = 1e12))$par
    }

    # Update prior distribution
    raw_prior <- t(est_clprior(min_missing, max_missing, itt, itt2, N2, num_treated, trt_friends, num_friends, p_prob, q_prob, alpha, beta))
    clprior_est <- normalize(raw_prior)

    if(abs(log_lik - old_log_lik) > epsilon) {
      old_log_lik = log_lik
    } else {
      results$means0 = sapply(1:4, function(i) sum((expit(params_est[1,i] + num_friends_x*params_est[2,i]) *
          (clpost_est[,1:xlen] + clpost_est[,1:xlen + xlen] + clpost_est[,1:xlen + 2*xlen] + clpost_est[,1:xlen + 3*xlen]))[num_friends_x >= 0])/
          sum((num_friends_x >= 0) * (clpost_est[,1:xlen] + clpost_est[,1:xlen + xlen] + clpost_est[,1:xlen + 2*xlen] + clpost_est[,1:xlen + 3*xlen])))
      results$params = params_est
      results$p_prob = p_prob
      results$q_prob = q_prob
      results$log_lik = log_lik
      results$iter = iter
      results$eff_obs = c(sum((num_friends_x >= 0)*clpost_est[,1:xlen]), sum((num_friends_x >= 0)*clpost_est[,xlen+1:xlen]),
                          sum((num_friends_x >= 0)*clpost_est[,2*xlen+1:xlen]), sum((num_friends_x >= 0)*clpost_est[,3*xlen+1:xlen]))
    }
  }

  return(results)
}










# Proposed EM estimator (logistic) with in- and out-of-village covariate
########################################################################

est_effects_em_binom_io <- function(data, p_prob_in, q_prob_in, alpha_in, beta_in, density_in,
                           p_prob_out, q_prob_out, alpha_out, beta_out, density_out,
                           min_missing = -10, max_missing = 10, max_iter = 200, epsilon = 0.001) {
  N_in <- data$N_in
  N_out <- data$N_out
  num_correct <- data$num_correct
  num_attempt <- data$num_attempt
  num_treated <- data$num_treated
  num_treated_in <- data$num_treated_in
  num_treated_out <- data$num_treated_out
  itt <- data$itt
  itt2 <- data$itt2
  trt_friends <- data$trt_friends
  trt_friends_in <- data$trt_friends_in
  trt_friends_out <- data$trt_friends_out
  num_friends <- data$num_friends
  num_friends_in <- data$num_friends_in
  num_friends_out <- data$num_friends_out

  results <- list()

  num_friends_x = sapply(min_missing:max_missing, function(x) num_friends + x)
  xlen = max_missing - min_missing + 1
  old_log_lik = -Inf

  # Estimate prior distribution of exposure conditions x true degrees (lexographical order)
  raw_prior <- t(est_clprior_io(min_missing, max_missing, itt, itt2, N_in, N_out, num_treated_in, num_treated_out, trt_friends_in, trt_friends_out,
                                num_friends_in, num_friends_out, p_prob_in, q_prob_in, p_prob_out, q_prob_out, alpha_in, beta_in, alpha_out, beta_out))
  clprior_est <- normalize(raw_prior)

  # Initialize parameters based on prior
  params_est <- matrix(data = 0, nrow = 2, ncol = 4)
  for(cl in 0:3) {
    params_est[,cl + 1] = optim(par = params_est[,cl+1], fn = param_binom_objective, post_mat = clprior_est[,1:xlen+cl*xlen],
                                num_friends_x = num_friends_x, n = num_attempt, a = num_correct, method = "BFGS", control = list(factr = 1e12))$par
  }

  for(iter in 1:max_iter) {
    # Estimate posterior distribution of exposure conditions x true degrees
    raw_post <- est_clpost_binom(clprior_est, num_friends_x, n = num_attempt, a = num_correct, params_est)
    log_lik <- sum(log(rowSums(raw_post)))
    clpost_est <- normalize(raw_post)

    # Marginal posterior of true degrees
    post_mat = t(clpost_est[,1:xlen] + clpost_est[,xlen + 1:xlen] + clpost_est[,2*xlen + 1:xlen] + clpost_est[,3*xlen + 1:xlen])

    # Estimate corruption in network p_in, p_out, q_in, and q_out
    opt = optim(fn = corr_objective_io, par = c(p_prob_in,q_prob_in, p_prob_out, q_prob_out), post_mat = post_mat,
                num_friends_in = num_friends_in, num_friends_out = num_friends_out, N_in = N_in, N_out = N_out,
                min_missing = min_missing, max_missing = max_missing, alpha_in = alpha_in, beta_in = beta_in, alpha_out = alpha_out, beta_out = beta_out,
                method = "L-BFGS-B", lower = c(1e-06,1e-06,1e-06,1e-10), upper = c(0.7,density_in*5,0.5,density_out*5), control = list(factr = 1e12))
    p_prob_in = opt$par[1]
    q_prob_in = opt$par[2]
    p_prob_out = opt$par[3]
    q_prob_out = opt$par[4]

    # Estimate parameters
    for(cl in 0:3) {
      params_est[,cl + 1] = optim(par = params_est[,cl+1], fn = param_binom_objective, post_mat = clpost_est[,1:xlen+cl*xlen],
                                  num_friends_x = num_friends_x, n = num_attempt, a = num_correct, method = "BFGS", control = list(factr = 1e12))$par
    }

    # Update prior distribution
    raw_prior <- t(est_clprior_io(min_missing, max_missing, itt, itt2, N_in, N_out, num_treated_in, num_treated_out, trt_friends_in, trt_friends_out,
                   num_friends_in, num_friends_out, p_prob_in, q_prob_in, p_prob_out, q_prob_out, alpha_in, beta_in, alpha_out, beta_out))
    clprior_est <- normalize(raw_prior)

    if(!is.na(log_lik) &  abs(log_lik - old_log_lik) > epsilon) {
      old_log_lik = log_lik
    } else {
      results$means0 = sapply(1:4, function(i) sum((expit(params_est[1,i] + num_friends_x*params_est[2,i]) *
          (clpost_est[,1:xlen] + clpost_est[,1:xlen + xlen] + clpost_est[,1:xlen + 2*xlen] + clpost_est[,1:xlen + 3*xlen]))[num_friends_x >= 0])/
          sum((num_friends_x >= 0) * (clpost_est[,1:xlen] + clpost_est[,1:xlen + xlen] + clpost_est[,1:xlen + 2*xlen] + clpost_est[,1:xlen + 3*xlen])))
      results$params = params_est
      results$p_prob_in = p_prob_in
      results$q_prob_in = q_prob_in
      results$p_prob_out = p_prob_out
      results$q_prob_out = q_prob_out
      results$log_lik = log_lik
      results$iter = iter
      results$eff_obs = c(sum((num_friends_x >= 0)*clpost_est[,1:xlen]), sum((num_friends_x >= 0)*clpost_est[,xlen+1:xlen]),
                          sum((num_friends_x >= 0)*clpost_est[,2*xlen+1:xlen]), sum((num_friends_x >= 0)*clpost_est[,3*xlen+1:xlen]))

    }
  }

  return(results)
}
