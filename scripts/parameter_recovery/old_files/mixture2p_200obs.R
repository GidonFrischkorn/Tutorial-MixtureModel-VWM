# Performs parameter recovery simulations for 2-parameter mixture model with fixed population
# parameters and 100 observations per subject
# The figures in https://github.com/GidonFrischkorn/Tutorial-MixtureModel-VWM/issues/16#issuecomment-2053756837 were generated using this script

library(here)
library(glue)
source(here("scripts/parameter_recovery/mixture2p_functions.R"))

# Parameters
N_replicates <- 250
N_subj <- 40
N_obs <- 200
range_thetat_mu <- qlogis(c(0.3,0.9))
thetat_sd <- 0.3
rnage_log_kappa <- log(c(2,15))
log_kappa_sd <- 0.2



i <- 0
fits <- replicate(N_replicates, {
  i <<- i + 1
  thetat_mu <- runif(1, range_thetat_mu[1], range_thetat_mu[2])
  log_kappa <- runif(1, rnage_log_kappa[1], rnage_log_kappa[2])
  print(paste0('Replicate ', i))
  fit_ml_and_bmm(N_subj, N_obs = N_obs, thetat_mu, thetat_sd, log_kappa, log_kappa_sd, cores = 4)
}, simplify = FALSE)
saveRDS(fits, here(glue('output/par_rec_fits{N_obs}.rds')))
