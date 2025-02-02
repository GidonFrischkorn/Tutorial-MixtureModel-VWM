library(here)
library(glue)
source(here("scripts/parameter_recovery/mixture2p_functions.R"))

# Parameters
set.seed(124)
N_replicates <- 200
N_subj <- 150
N_obs <- 50
thetat_mu <- 0
thetat_sd <- 0.4
log_kappa <- 1
log_kappa_sd <- 0.2



i <- 0
fits <- replicate(N_replicates, {
  i <<- i + 1
  print(paste0('Replicate ', i))
  fit_ml_and_bmm(N_subj, N_obs = N_obs, thetat_mu, thetat_sd, log_kappa, log_kappa_sd, cores = 4)
}, simplify = FALSE)
saveRDS(fits, here(glue('output/par_rec_fits{N_obs}_{N_subj}subj_{thetat_mu}theta_{log_kappa}log_kappa.rds')))
