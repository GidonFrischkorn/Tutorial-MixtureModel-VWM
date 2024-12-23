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
N_replicates <- 200
cond_N_subj <- c(20,40,80)
cond_N_obs <- c(25,50,100,200)
range_thetat_mu <- qlogis(c(0.3,0.95))
thetat_sd <- 0.3
rnage_log_kappa <- log(c(2,15))
log_kappa_sd <- 0.2

df_conditions <- expand.grid(
  n_sub = cond_N_subj,
  n_trials = cond_N_obs
)

## Report on progress automatically
handlers("cli")

## Parallelize on local machine
plan(multisession)

## See https://progressr.futureverse.org/#a-progressor-cannot-be-created-in-the-global-environment
## for why we use local() here
reprex <- local({
  p <- progressor(steps = nrow(df_conditions) * N_replicates)
  y <- foreach(i = 1:N_replicates, .combine='c', .multicombine=TRUE) %dofuture% {
    
    for(j in 1:nrow(df_conditions)){
      N_sub <- df_conditions$n_sub[j]
      N_obs <- df_conditions$n_trials[j]
      
      thetat_mu <- runif(1, range_thetat_mu[1], range_thetat_mu[2])
      log_kappa <- runif(1, rnage_log_kappa[1], rnage_log_kappa[2])
      results_list <- fit_ml_and_bmm(N_sub, N_obs = N_obs, thetat_mu, thetat_sd, log_kappa, log_kappa_sd, cores = 1)
      
      saveRDS(results_list, here("output","recovery_results",glue('par_rec_fits_{N_sub}sub_{N_obs}trials_rep{i}.rds')))
      
      ## Report on progress thus far with a custom message
      p(sprintf("(j, i) = (%d, %d)", j, i))
    }
  }
})

