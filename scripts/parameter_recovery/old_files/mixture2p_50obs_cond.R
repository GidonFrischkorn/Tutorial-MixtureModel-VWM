library(here)
library(glue)
source(here("scripts/parameter_recovery/mixture2p_functions.R"))

# Parameters
set.seed(124)
N_replicates <- 5
N_subj <- 40
N_obs <- 50
thetat_mu <- 1
thetat_sd <- 0.4
log_kappa <- 3
log_kappa_sd <- 0.2
log_kappa_cond_mu <- -0.2
log_kappa_cond_sd <- 0.00000001


i <- 0
fits <- replicate(N_replicates, {
  i <<- i + 1
  print(paste0('Replicate ', i))
  fit_ml_and_bmm_cond_diff(N_subj, N_obs = N_obs, thetat_mu, thetat_sd, log_kappa, log_kappa_sd, 
                           log_kappa_cond_mu = log_kappa_cond_mu,
                           log_kappa_cond_sd = log_kappa_cond_sd)
}, simplify = FALSE)
saveRDS(fits, here(glue('output/par_rec_fits{N_obs}.rds')))





# Parameters
set.seed(124)
N_replicates <- 5
N_subj <- 40
N_obs <- 50
thetat_mu <- 1
thetat_sd <- 0.4
thetat_cond_mu <- 0.5
thetat_cond_sd <- 0.3
log_kappa <- 3
log_kappa_sd <- 0.2


i <- 0
fits_theta <- replicate(N_replicates, {
  i <<- i + 1
  print(paste0('Replicate ', i))
  fit_ml_and_bmm_cond_diff(N_subj, N_obs = N_obs, thetat_mu, thetat_sd, log_kappa, log_kappa_sd, 
                           thetat_cond_mu = thetat_cond_mu,
                           thetat_cond_sd = thetat_cond_sd)
}, simplify = FALSE)