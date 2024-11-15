# Performs parameter recovery simulations for 2-parameter mixture model with fixed population
# parameters and 500 observations per subject
# The figures in https://github.com/GidonFrischkorn/Tutorial-MixtureModel-VWM/issues/16#issuecomment-2053756837 were generated using this script

library(here)
library(glue)
source(here("scripts/parameter_recovery/mixture2p_functions.R"))

# Parameters
set.seed(124)
N_replicates <- 200
N_subj <- 40
N_obs <- 500
thetat_mu <- 1
thetat_sd <- 0.4
log_kappa <- 3
log_kappa_sd <- 0.2



i <- 0
fits <- replicate(N_replicates, {
  i <<- i + 1
  print(paste0('Replicate ', i))
  fit_ml_and_bmm(N_subj, N_obs = N_obs, thetat_mu, thetat_sd, log_kappa, log_kappa_sd)
}, simplify = FALSE)
saveRDS(fits, here(glue('output/par_rec_fits{N_obs}.rds')))
