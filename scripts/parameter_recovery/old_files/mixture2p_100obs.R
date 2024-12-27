# Performs parameter recovery simulations for 2-parameter mixture model with fixed population
# parameters and 100 observations per subject
# The figures in https://github.com/GidonFrischkorn/Tutorial-MixtureModel-VWM/issues/16#issuecomment-2053756837 were generated using this script

library(here)
library(glue)
library(foreach)
library(doFuture)
library(progressr)
source(here("scripts/parameter_recovery/mixture2p_functions.R"))

# Parameters
N_replicates <- 10
N_subj <- 40
cond_N_obs <- c(25,50,100,200)
range_thetat_mu <- qlogis(c(0.3,0.9))
thetat_sd <- 0.3
rnage_log_kappa <- log(c(2,15))
log_kappa_sd <- 0.2

# setup parallel backend to use many processors
cores = detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

## Report on progress automatically
handlers(global = TRUE)

## Parallelize on local machine
plan(multisession)

foreach::foreach (j = 1:length(cond_N_obs)) %dopar% {
  results_list <- list()
  N_obs <- cond_N_obs[j]
  
  for (i in 1:N_replicates){
    thetat_mu <- runif(1, range_thetat_mu[1], range_thetat_mu[2])
    log_kappa <- runif(1, rnage_log_kappa[1], rnage_log_kappa[2])
    results_list[[i]] <- fit_ml_and_bmm(N_subj, N_obs = N_obs, thetat_mu, thetat_sd, log_kappa, log_kappa_sd, cores = 1)
  }
  
  saveRDS(results_list, here(glue('output/par_rec_fits_{N_sub}sub_{N_obs}trials.rds')))
}

#stop cluster
stopCluster(cl)
