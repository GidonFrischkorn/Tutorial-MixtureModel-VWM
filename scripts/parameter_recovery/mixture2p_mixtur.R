library(pacman)
p_load(here)
p_load_current_gh(
  "venpopov/mixtur@rescale_parameters", "venpopov/bmm",
  "djnavarro/queue"
)

# number of cores to use for the simulation
cores <- 30
if (cores > future::availableCores()) cores <- future::availableCores() - 1

# generates data from the mixture2p model, and fits it via the mixtur package
fit_m2_ml <- function(n_obs, pars) {
  pars$id <- seq_len(nrow(pars))
  pars <- pars[, c("id", "pmem", "kappa")]
  N_subj <- nrow(pars) # nolint: object_name_linter.

  data <- list()
  for (i in 1:N_subj) {
    data[[i]] <- bmm::rmixture2p(
      n_obs,
      mu = 0,
      pMem = pars$pmem[i],
      kappa = pars$kappa[i]
    )
  }

  data <- data.frame(
    id = rep(1:N_subj, each = n_obs),
    y = unlist(data),
    target = rep(0, N_subj * n_obs)
  )

  par_ml <- mixtur::fit_mixtur(
    data,
    model = "2_component",
    unit = "radians",
    response_var = "y"
  )

  names(par_ml) <- c("id", "kappa", "pmem", "pguess")
  par_ml <- par_ml[, c("id", "pmem", "kappa")]
  par_ml$Nobs <- n_obs
  par_ml
}

# store results
file <- here("output/mixture2p_mixtur_fullgrid_recovery_unconstrained.RData")

if (!file.exists(file)) {
  # setup simulation details
  pars <- expand.grid(
    kappa = seq(1, 16, 0.5),
    pmem = seq(0.1, 1, 0.025)
  )
  n_obs <- c(20, 50, 100, 200, 500) # number of observations per design cell
  rep <- 1000 # number of replications of each simulation

  # setup parallel tasks in external r session via the queue package
  queue <- Queue$new(workers = cores)
  for (N in n_obs) {
    for (r in 1:rep) {
      queue$add(fit_m2_ml, args = list(n_obs = N, pars = pars))
    }
  }

  out <- queue$run(message = "verbose")
  res <- out$result

  dir <- dirname(file)
  if (!dir.exists(dir)) dir.create(dir)
  save(res, pars, file = file)
}
