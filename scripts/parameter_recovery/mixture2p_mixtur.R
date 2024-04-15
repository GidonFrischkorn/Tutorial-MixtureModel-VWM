library(bmm)
library(mixtur)
library(here)

fit_m2_ml <- function(N_obs) {
  pars <- expand.grid(kappa = seq(1, 16, 0.5),
                      pmem = seq(0.6, 1, 0.025))
  pars$id <- 1:nrow(pars)
  pars <- pars[,c('id','pmem','kappa')]
  N_subj <- nrow(pars)
  
  data <- list()
  for (i in 1:N_subj) {
    data[[i]] <- rmixture2p(N_obs, mu = 0, pMem = pars$pmem[i], kappa = pars$kappa[i])
  }
  
  data <- data.frame(id = rep(1:N_subj, each = N_obs), 
                     y = unlist(data),
                     target = rep(0, N_subj * N_obs))
  
  par_ml <- fit_mixtur(data, model = "2_component", unit = "radians", response_var = "y")
  names(par_ml) <- c('id','kappa','pmem','pguess')
  par_ml <- par_ml[,c('id','pmem','kappa')]
  par_ml
}


N_obs <- c(20, 50, 200, 500)

fits20 <- replicate(200, fit_m2_ml(N_obs = 20))
fits50 <- replicate(200, fit_m2_ml(N_obs = 50))
fits100 <- replicate(200, fit_m2_ml(N_obs = 100))
fits200 <- replicate(200, fit_m2_ml(N_obs = 200))
fits500 <- replicate(200, fit_m2_ml(N_obs = 500))

save(fits20, fits50, fits100, fits200, fits500, file = here('output/mixture2p_mixtur_grid_recovery.RData'))
