# THIS SCRIPT NEVER WORKED, AS IT WAS UNFINISHED 

# library(pacman)
# p_load(here)
# p_load_current_gh("venpopov/mixtur@efficiency", "venpopov/bmm", 'djnavarro/queue')

# # number of cores to use for the simulation
# cores = 28
# if (cores > future::availableCores()) cores = future::availableCores()-1

# # generates data from the mixture2p model, and fits it via the mixtur package
# fit_m3_ml <- function(N_obs, pars) {
#   pars$id <- 1:nrow(pars)
#   pars <- pars[, c('id', 'pmem', 'pnt', 'kappa')]
#   N_subj <- nrow(pars)
  
#   data <- list()
#   for (i in 1:N_subj) {
#     data[[i]] <- bmm::rmixture3p(
#       N_obs,
#       mu = 0,
#       pMem = pars$pmem[i],
#       kappa = pars$kappa[i],
#       pNT = pars$pnt[i]
#     )
#   }
  
#   data <- data.frame(
#     id = rep(1:N_subj, each = N_obs),
#     y = unlist(data),
#     target = rep(0, N_subj * N_obs)
#   )
  
#   par_ml <- mixtur::fit_mixtur(
#     data,
#     model = "2_component",
#     unit = "radians",
#     response_var = "y"
#   )
  
#   names(par_ml) <- c('id', 'kappa', 'pmem', 'pguess')
#   par_ml <- par_ml[, c('id', 'pmem', 'kappa')]
#   par_ml$Nobs <- N_obs
#   par_ml
# }

# # store results
# file = here('output/mixture3p_mixtur_fullgrid_recovery.RData')

# if (!file.exists(file)) {
#   # setup simulation details
#   pars <- expand.grid(
#     kappa = seq(1, 16, 0.5), 
#     pmem = seq(0.1, 1, 0.05),
#     pnt = seq(0, 1, 0.05)
#   ) 
#   pars <- dplyr::filter(pars, pmem + pnt <= 1)
#   N_obs <- c(20, 50, 100, 200, 500)       # number of observations per design cell
#   rep <- 1000                             # number of replications of each simulation
  
#   # setup parallel tasks in external r session via the queue package
#   queue <- Queue$new(workers = cores)
#   for (N in N_obs) {
#     for (r in 1:rep) {
#       queue$add(fit_m3_ml, args = list(N_obs = N, pars = pars))  
#     }
#   } 
  
#   out <- queue$run(message = "verbose")
#   res <- out$result
  
#   dir <- dirname(file)
#   if (!dir.exists(dir)) dir.create(dir)
#   save(res, pars, file = file)
# }



