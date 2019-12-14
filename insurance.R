# Setup
#######

library(foreign)

setwd("insurance/data")

raw_net_data <- read.dta("0422allinforawnet.dta")
survey_data <- read.dta("0422survey.dta")

expit = function(x) exp(x)/(1+exp(x))



# Create network data (mirroring insurance/do/rawnet.do)
########################################################

N2 = 4662
villages = aggregate(raw_net_data$village, by=list(raw_net_data$id), FUN = function(x) x[1])

raw_net_data = raw_net_data[raw_net_data$network_missname != 1 | !is.na(raw_net_data$network_id),]
raw_net_data$network_pre_intensive = (raw_net_data$delay == 0) & (raw_net_data$intensive == 1)
raw_net_data$mutual = sapply(1:dim(raw_net_data)[1],
                             function(x) any(raw_net_data$id == raw_net_data$network_id[x] & raw_net_data$network_id == raw_net_data$id[x])*1)

# "Normal" network
network_data = aggregate(rep(1,length(raw_net_data$id)), by=list(raw_net_data$id), FUN = sum, na.rm = TRUE)
network_data_in = aggregate(raw_net_data$village == raw_net_data$network_village, by=list(raw_net_data$id), FUN = sum, na.rm = TRUE)
network_data_out = aggregate(raw_net_data$village != raw_net_data$network_village, by=list(raw_net_data$id), FUN = sum, na.rm = TRUE)

num_trt_friends = aggregate(raw_net_data$network_pre_intensive, by=list(raw_net_data$id), FUN = sum, na.rm = TRUE)
num_trt_friends_in = aggregate((raw_net_data$village == raw_net_data$network_village) & raw_net_data$network_pre_intensive,
                               by=list(raw_net_data$id), FUN = sum, na.rm = TRUE)
num_trt_friends_out = aggregate((raw_net_data$village != raw_net_data$network_village) & raw_net_data$network_pre_intensive,
                                by=list(raw_net_data$id), FUN = sum, na.rm = TRUE)

# Network of mutual friends
num_mutual_friends = aggregate(raw_net_data$mutual, by=list(raw_net_data$id), FUN = sum, na.rm = TRUE)
num_mutual_friends_in = aggregate(raw_net_data$mutual & (raw_net_data$village == raw_net_data$network_village),
                                  by=list(raw_net_data$id), FUN = sum, na.rm = TRUE)
num_mutual_friends_out = aggregate(raw_net_data$mutual & (raw_net_data$village != raw_net_data$network_village),
                                   by=list(raw_net_data$id), FUN = sum, na.rm = TRUE)

num_mutual_trt_friends = aggregate(raw_net_data$network_pre_intensive*raw_net_data$mutual, by=list(raw_net_data$id), FUN = sum, na.rm = TRUE)
num_mutual_trt_friends_in = aggregate(raw_net_data$network_pre_intensive*raw_net_data$mutual & (raw_net_data$village == raw_net_data$network_village),
                                      by=list(raw_net_data$id), FUN = sum, na.rm = TRUE)
num_mutual_trt_friends_out = aggregate(raw_net_data$network_pre_intensive*raw_net_data$mutual & (raw_net_data$village != raw_net_data$network_village),
                                       by=list(raw_net_data$id), FUN = sum, na.rm = TRUE)

# Network of weak ties
temp_net_data1 = raw_net_data[,c("id","network_id","village")]
temp_net_data2 = raw_net_data[,c("id","network_id","network_pre_intensive","network_village")]
temp_net_data = merge(temp_net_data1, temp_net_data2, by.x = "network_id", by.y = "id")
temp_net_data$network_id = NULL
colnames(temp_net_data)[3] = "network_id"
temp_net_data1 = rbind(temp_net_data, raw_net_data[,c("id","village","network_id","network_pre_intensive","network_village")])

temp_net_data = temp_net_data1[!duplicated(temp_net_data1[,c(1:3)]),]
num_weak_friends = aggregate(rep(1,length(temp_net_data$id)), by=list(temp_net_data$id), FUN = sum, na.rm = TRUE)
num_weak_friends_in = aggregate(temp_net_data$village == temp_net_data$network_village, by=list(temp_net_data$id), FUN = sum, na.rm = TRUE)
num_weak_friends_out = aggregate(temp_net_data$village != temp_net_data$network_village, by=list(temp_net_data$id), FUN = sum, na.rm = TRUE)

num_weak_trt_friends = aggregate(temp_net_data$network_pre_intensive, by=list(temp_net_data$id), FUN = sum, na.rm = TRUE)
num_weak_trt_friends_in = aggregate(temp_net_data$network_pre_intensive & (temp_net_data$village == temp_net_data$network_village),
                                    by=list(temp_net_data$id), FUN = sum, na.rm = TRUE)
num_weak_trt_friends_out = aggregate(temp_net_data$network_pre_intensive & (temp_net_data$village != temp_net_data$network_village),
                                     by=list(temp_net_data$id), FUN = sum, na.rm = TRUE)



#merge network data altogther
colnames(network_data)[2] = "num_friends"
colnames(network_data)[1] = "id"

network_data$num_friends_in = network_data_in$x
network_data$num_friends_out = network_data_out$x
network_data$num_trt_friends = num_trt_friends$x
network_data$num_trt_friends_in = num_trt_friends_in$x
network_data$num_trt_friends_out = num_trt_friends_out$x

network_data$num_mutual_friends = num_mutual_friends$x
network_data$num_mutual_friends_in = num_mutual_friends_in$x
network_data$num_mutual_friends_out = num_mutual_friends_out$x
network_data$num_mutual_trt_friends = num_mutual_trt_friends$x
network_data$num_mutual_trt_friends_in = num_mutual_trt_friends_in$x
network_data$num_mutual_trt_friends_out = num_mutual_trt_friends_out$x

network_data$num_weak_friends = num_weak_friends$x
network_data$num_weak_friends_in = num_weak_friends_in$x
network_data$num_weak_friends_out = num_weak_friends_out$x
network_data$num_weak_trt_friends = num_weak_trt_friends$x
network_data$num_weak_trt_friends_in = num_weak_trt_friends_in$x
network_data$num_weak_trt_friends_out = num_weak_trt_friends_out$x

network_data$village = villages$x
network_data$trt = sapply(1:N2, function(i) any(raw_net_data$network_pre_intensive[raw_net_data$network_id == network_data$id[i]],na.rm = TRUE))
for(i in 1:dim(survey_data)[1]) {
  if(survey_data$id[i] %in% network_data$id) {
    network_data$trt[which(network_data$id == survey_data$id[i])] = survey_data$intensive[i] == 1 & survey_data$delay[i] == 0
  }
}
network_data$N_in = sapply(1:N2, function(i) sum(network_data$village[i] == network_data$village))
network_data$N_out = dim(network_data)[1] - network_data$N_in
network_data$num_treated_in = sapply(1:N2, function(i) sum((network_data$village[i] == network_data$village) & network_data$trt))
network_data$num_treated_out = sum(network_data$trt) - network_data$num_treated_in



#Merge datasets
data <- merge(network_data, survey_data, by ="id")










# Estimation
############

# Restrict to subjects in the second round of sessions
data = data[data$delay == 1 & data$info_none == 1,]
N = dim(data)[1]

insurance_data <- list()
insurance_data$N = N
insurance_data$N2 = N2
insurance_data$num_correct = round(data$understanding*5)
insurance_data$num_attempt = rep(5,length(insurance_data$num_correct))
insurance_data$num_treated = rep(1003,N)
insurance_data$itt = data$intensive & data$delay == 1 #Subjects who attended an intensive session in round 2
insurance_data$itt2 = data$intensive & data$delay == 0 #Subjects who attended an intensive session in round 1
insurance_data$trt_friends = data$num_trt_friends
insurance_data$num_friends = data$num_friends
#insurance_data$trt_friends = data$num_mutual_trt_friends
#insurance_data$num_friends = data$num_mutual_friends
#insurance_data$trt_friends = data$num_weak_trt_friends
#insurance_data$num_friends = data$num_weak_friends



density = mean(data$num_friends)/(N2-1)
rho_est = (var(data$num_friends)/(N2-1)/density/(1-density)-1)/(N-1)
alpha = density*(1/rho_est-1)
beta = (1-density)*(1/rho_est-1)

init_values = expand.grid(p_prob = c(0.05,0.1, 0.15), q_prob = c(0.2,0.3,0.4)*density)

results <- list()

for(i in 1:dim(init_values)[1]) {
  p_prob = init_values[i,1]
  q_prob = init_values[i,2]

  results[[i]] = est_effects_em_binom(insurance_data, p_prob, q_prob, alpha, beta, min_missing = -10, max_missing = 5, max_iter = 200, epsilon = 0.001)
  results[[i]]$init_p_prob = p_prob
  results[[i]]$init_q_prob = q_prob

}

save(results,file="insurance_results.RData")





# Bootstrap standard errors
index = 5

min_missing = -10
max_missing = 5

p_prob = results[[index]]$init_p_prob
q_prob = results[[index]]$init_q_prob

num_treated <- insurance_data$num_treated
itt <- insurance_data$itt
itt2 <- insurance_data$itt2
trt_friends <- insurance_data$trt_friends
num_friends <- insurance_data$num_friends

num_friends_x = sapply(min_missing:max_missing, function(x) num_friends + x)
xlen = max_missing - min_missing + 1
raw_prior <- t(est_clprior(min_missing, max_missing, itt, itt2, N2, num_treated, trt_friends, num_friends, p_prob, q_prob, alpha, beta))
clprior_est <- normalize(raw_prior)
clprior_cs = t(apply(clprior_est, 1, cumsum))

boot_results = list()
for(i in 1:50) {
  # Draw true degree and exposure condition
  U <- runif(N)
  latent_var = (sapply(1:N, function(i) which.min(clprior_cs[i,] < U[i])) - 1)
  num_friends_latent = transform(num_friends + (latent_var  %% (max_missing - min_missing + 1)) + min_missing)
  exposure = floor(latent_var/(max_missing - min_missing + 1)) + 1

  insurance_data$num_correct <- rbinom(N, 5, expit(results[[index]]$params[1,exposure] + results[[index]]$params[2,exposure]*num_friends_latent))

  boot_results[[i]] = est_effects_em_binom(insurance_data, p_prob, q_prob, alpha, beta, min_missing = -10, max_missing = 5, max_iter = 200, epsilon = 0.001)
}

save(boot_results,file="insurance_boot_results.RData")










## With village covariates
##########################

insurance_data$num_correct = round(data$understanding*5)
insurance_data$N_in = data$N_in
insurance_data$N_out = data$N_out

insurance_data$num_treated_in = data$num_treated_in
insurance_data$num_treated_out = data$num_treated_out

insurance_data$trt_friends_in = data$num_trt_friends_in
insurance_data$trt_friends_out = data$num_trt_friends_out
insurance_data$num_friends_in = data$num_friends_in
insurance_data$num_friends_out = data$num_friends_out
#insurance_data$trt_friends_in = data$num_mutual_trt_friends_in
#insurance_data$trt_friends_out = data$num_mutual_trt_friends_out
#insurance_data$num_friends_in = data$num_mutual_friends_in
#insurance_data$num_friends_out = data$num_mutual_friends_out
#insurance_data$trt_friends_in = data$num_weak_trt_friends_in
#insurance_data$trt_friends_out = data$num_weak_trt_friends_out
#insurance_data$num_friends_in = data$num_weak_friends_in
#insurance_data$num_friends_out = data$num_weak_friends_out



density_in = sum(data$num_friends_in)/sum(data$N_in)
density_out = sum(data$num_friends_out)/sum(data$N_out)
rho_est = (var(data$num_friends)/(N2-1)/density/(1-density)-1)/(N-1)
alpha_in = density_in*(1/rho_est-1)
beta_in = (1-density_in)*(1/rho_est-1)
alpha_out = density_out*(1/rho_est-1)
beta_out = (1-density_out)*(1/rho_est-1)

p_prob_in = 0.1
q_prob_in = density_in * 0.3
p_prob_out = 0.1
q_prob_out = density_out * 0.3

results = est_effects_em_binom_io(insurance_data, p_prob_in, q_prob_in, alpha_in, beta_in, density_in,
                                  p_prob_out, q_prob_out, alpha_out, beta_out, density_out,
                                  min_missing = -5, max_missing = 5, max_iter = 200, epsilon = 0.001)
results$init_p_prob_in = p_prob_in
results$init_q_prob_in = q_prob_in
results$init_p_prob_out = p_prob_out
results$init_q_prob_out = q_prob_out

save(results,file="insurance_results_io.RData")





# Bootstrap standard errors
min_missing = -5
max_missing = 5

N_in <- insurance_data$N_in
N_out <- insurance_data$N_out
num_treated_in <- insurance_data$num_treated_in
num_treated_out <- insurance_data$num_treated_out
itt <- insurance_data$itt
itt2 <- insurance_data$itt2
trt_friends_in <- insurance_data$trt_friends_in
trt_friends_out <- insurance_data$trt_friends_out
num_friends <- insurance_data$num_friends
num_friends_in <- insurance_data$num_friends_in
num_friends_out <- insurance_data$num_friends_out

num_friends_x = sapply(min_missing:max_missing, function(x) num_friends + x)
xlen = max_missing - min_missing + 1
raw_prior <- t(est_clprior_io(min_missing, max_missing, itt, itt2, N_in, N_out, num_treated_in, num_treated_out, trt_friends_in, trt_friends_out,
                              num_friends_in, num_friends_out, p_prob_in, q_prob_in, p_prob_out, q_prob_out, alpha_in, beta_in, alpha_out, beta_out))
clprior_est <- normalize(raw_prior)
clprior_cs = t(apply(clprior_est, 1, cumsum))

boot_results = list()
for(i in 1:50) {
  # Draw true degree and exposure condition
  U <- runif(N)
  latent_var = (sapply(1:N, function(i) which.min(clprior_cs[i,] < U[i])) - 1)
  num_friends_latent = transform(num_friends + (latent_var  %% (max_missing - min_missing + 1)) + min_missing)
  exposure = floor(latent_var/(max_missing - min_missing + 1)) + 1

  insurance_data$num_correct <- rbinom(N, 5, expit(results[[index]]$params[1,exposure] + results[[index]]$params[2,exposure]*num_friends_latent))

  boot_results[[i]] = est_effects_em_binom_io(insurance_data, p_prob_in, q_prob_in, alpha_in, beta_in, density_in,
                                  p_prob_out, q_prob_out, alpha_out, beta_out, density_out,
                                  min_missing = -5, max_missing = 5, max_iter = 200, epsilon = 0.001)
}

save(boot_results,file="insurance_boot_results_io.RData")










# Basic GLM model assuming no mismeasurement
############################################

exposure = as.factor(itt*2 + (trt_friends>0))
model <- glm(y ~ exposure + exposure*num_friends - 1 - num_friends, family = binomial(link = "logit"), weights = rep(5,length(y)))
summary(model)

means0 = rep(0,4)
means0[1] = mean(expit(model$coefficients[1] + model$coefficients[5]*num_friends))
means0[2] = mean(expit(model$coefficients[2] + model$coefficients[6]*num_friends))
means0[3] = mean(expit(model$coefficients[3] + model$coefficients[7]*num_friends))
means0[4] = mean(expit(model$coefficients[4] + model$coefficients[8]*num_friends))

boot_glm <- matrix(data = 0, nrow = 100, ncol = 4)
for(i in 1:100) {
  new_y = rbinom(length(y),5,model$fitted.values)
  new_model <- glm(new_y/5 ~ exposure + exposure*num_friends - 1 - num_friends, family = binomial(link = "logit"), weights = rep(5,length(y)))
  boot_glm[i,1] = mean(expit(new_model$coefficients[1] + new_model$coefficients[5]*num_friends))
  boot_glm[i,2] = mean(expit(new_model$coefficients[2] + new_model$coefficients[6]*num_friends))
  boot_glm[i,3] = mean(expit(new_model$coefficients[3] + new_model$coefficients[7]*num_friends))
  boot_glm[i,4] = mean(expit(new_model$coefficients[4] + new_model$coefficients[8]*num_friends))
}

save(boot_glm, means, model, file = "insurance_results_comp.RData")

means[3] - means[1] # Effect of intensive session (direct effect)
means[2] - means[1] # Network intensive (indirect effect)
means[4] + means[1] - means[2] - means[3] # Interaction

sd(sapply(1:100, function(i) boot_glm[i,3] - boot_glm[i,1]))
sd(sapply(1:100, function(i) boot_glm[i,2] - boot_glm[i,1]))
sd(sapply(1:100, function(i) boot_glm[i,4] + boot_glm[i,1] - boot_glm[i,2] - boot_glm[i,3]))
