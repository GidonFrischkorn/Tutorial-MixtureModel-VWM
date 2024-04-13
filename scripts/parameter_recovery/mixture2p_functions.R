library(mixtur)
library(bmm)

fit_ml_and_bmm <- function(N_subj, N_obs, thetat_mu, thetat_sd, log_kappa, log_kappa_sd, cores = 4) {
  thetats = rnorm(N_subj, thetat_mu, thetat_sd)
  log_kappas = rnorm(N_subj, log_kappa, log_kappa_sd)
  pars <- data.frame(id = 1:N_subj, pmem = plogis(thetats), kappa = exp(log_kappas))
  
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
  
  fit_bmm <- bmm(bmf(thetat ~ 1 + (1|id),
                     kappa ~ 1 + (1|id)),
                 data = data,
                 model = mixture2p('y'),
                 cores = cores,
                 backend = 'cmdstanr',
                 refresh = 0)
  
  fe_bmm <- fixef(fit_bmm)[c('kappa_Intercept','thetat_Intercept'),1]
  re_bmm <- ranef(fit_bmm)$id[,1,]
  par_bmm <- as.data.frame(t(t(re_bmm) + fe_bmm))
  par_bmm$id <- 1:N_subj
  names(par_bmm) <- c('kappa','pmem','id')
  par_bmm <- par_bmm[,c('id','pmem','kappa')]
  par_bmm$kappa <- exp(par_bmm$kappa)
  par_bmm$pmem <- exp(par_bmm$pmem)/(1+exp(par_bmm$pmem))
  
  return(list(pars = pars, par_ml = par_ml, par_bmm = par_bmm, fit_bmm = fit_bmm, data = data))
}

plot_par_rec <- function(fits) {
  with(fits, {
    par(mfrow = c(2,2))
    plot(pars$pmem, par_ml$pmem, xlab = 'True', ylab = 'Estimated (ML)', main = 'pmem')
    abline(0,1)
    plot(pars$kappa, par_ml$kappa, xlab = 'True', ylab = 'Estimated (ML)', main = 'kappa')
    abline(0,1)
    
    plot(pars$pmem, par_bmm$pmem, xlab = 'True', ylab = 'Estimated (BMM)', main = 'pmem')
    abline(0,1)
    plot(pars$kappa, par_bmm$kappa, xlab = 'True', ylab = 'Estimated (BMM)', main = 'kappa')
    abline(0,1)
  })
}


par_rec_stats <- function(fits) {
  with(fits, {
    data.frame(
      engine = c('ML','ML','BMM','BMM'),
      param = c('pmem','kappa','pmem','kappa'),
      r = c(cor.test(pars$pmem, par_ml$pmem)$estimate,
            cor.test(pars$kappa, par_ml$kappa)$estimate,
            cor.test(pars$pmem, par_bmm$pmem)$estimate,
            cor.test(pars$kappa, par_bmm$kappa)$estimate),
      rmse = c(sqrt(mean((pars$pmem - par_ml$pmem)^2)),
               sqrt(mean((pars$kappa - par_ml$kappa)^2)),
               sqrt(mean((pars$pmem - par_bmm$pmem)^2)),
               sqrt(mean((pars$kappa - par_bmm$kappa)^2))),
      pop_par_diff = c(mean(pars$pmem) - mean(par_ml$pmem),
                      mean(pars$kappa) - mean(par_ml$kappa),
                      mean(pars$pmem) - mean(par_bmm$pmem),
                      mean(pars$kappa) - mean(par_bmm$kappa))
    )
  })
}
