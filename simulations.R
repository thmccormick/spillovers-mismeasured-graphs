# Setup
#######

library(igraph)
library(ggplot2)
library(doSNOW)
library(foreach)
library(Rcpp)



# Load network data from Banerjee et al (2013)
village_indices = c(1:77)[-c(13,22)]
W <- list()
for(i in 1:75) {
  W1 <- as.matrix(read.csv(paste("villages/adj_borrowmoney_HH_vilno_",village_indices[i],".csv",sep =""), header = FALSE))
  W2 <- as.matrix(read.csv(paste("villages/adj_lendmoney_HH_vilno_",village_indices[i],".csv",sep =""), header = FALSE))
  W[[i]] = (W1 + W2) > 0
}



# Generate data for simulations
params_true <- list()
params_true$N = 250
params_true$base_prob = 0.25
params_true$mu = c(0,0.25,0.5,1)
params_true$deg_beta = c(0.05,0.1,0.05,0.1)
params_true$sigma2 = c(0.5,0.5,0.5,0.5)/2


source("methods.R")





# Simulations
#############

ps = seq(0,1/2,1/8)

for(q_name in 0:4) {
  q = q_name/8

  for(run in 1:10) {
    result <- list()

    for(i in 1:75) {
      # Simulate data
      params_true$N = dim(W[[i]])[1]
      degrees = rowSums(W[[i]])
      density = sum(W[[i]])/params_true$N/(params_true$N-1)
      q_prob = density*q
      sim_data <- gen_data_village(params_true, W[[i]], ps, q_prob)

      result[[i]] <- list()
      result[[i]]$truth = params_true$mu + mean(degrees[degrees > 0])*params_true$deg_beta
      result[[i]]$truth0 = params_true$mu + mean(degrees[degrees >= 0])*params_true$deg_beta

      # Estimation
      result[[i]]$em <- list()
      result[[i]]$as = lapply(sim_data, est_effects_ht, base_prob = params_true$base_prob)

      for(j in 1:length(ps)) {
        p_prob = ps[j]

        # Degree distribution parameters
        density_est = sum(sim_data[[j]]$W_0)/params_true$N/(params_true$N-1)
        rho_est = (var(rowSums(sim_data[[j]]$W_0))/(params_true$N-1)/density_est/(1-density_est)-1)/(params_true$N-2)
        alpha_est = density_est*(1/rho_est-1)
        beta_est = (1-density_est)*(1/rho_est-1)

        est_results = est_effects_em(sim_data[[j]], p_prob, q_prob, alpha_est, beta_est,
                        min_missing = -5+min(degrees - rowSums(sim_data[[j]]$W_0)), max_missing = 5+max(degrees - rowSums(sim_data[[j]]$W_0)))
        result[[i]]$em[[j]] = est_results
      }

      result[[i]]$time = stop - start
    }

    save(result, file = paste("simulation_results_q",q_name,"_",run,".RData",sep=""))
  }
}
