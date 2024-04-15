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

fit_ml_and_bmm_cond_diff <- function(N_subj, N_obs, thetat_mu, thetat_sd, log_kappa, log_kappa_sd, 
                                     thetat_cond_mu, thetat_cond_sd, log_kappa_cond_mu, 
                                     log_kappa_cond_sd, cores = 4) {
  thetats_A <- rnorm(N_subj, thetat_mu, thetat_sd)
  if (!missing(thetat_cond_mu)) {
    thetats_delta <- rnorm(N_subj, thetat_cond_mu, thetat_cond_sd)
  } else {
    thetats_delta <- rep(0, N_subj)
  }
  thetats_B <- thetats_A + thetats_delta
  
  log_kappas_A <- rnorm(N_subj, log_kappa, log_kappa_sd)
  if (!missing(log_kappa_cond_mu)) {
    log_kappas_delta <- rnorm(N_subj, log_kappa_cond_mu, log_kappa_cond_sd)
  } else {
    log_kappas_delta <- rep(0, N_subj)
  }
  log_kappas_B <- log_kappas_A + log_kappas_delta
  
  pars_A <- data.frame(id = 1:N_subj, pmem = plogis(thetats_A), kappa = exp(log_kappas_A), cond = "A")
  pars_B <- data.frame(id = 1:N_subj, pmem = plogis(thetats_B), kappa = exp(log_kappas_B), cond = "B")
  pars <- rbind(pars_A, pars_B)
  
  data <- list()
  for (i in 1:N_subj) {
    A <- rmixture2p(N_obs, mu = 0, pMem = pars[pars$cond == "A", ]$pmem[i], 
                    kappa = pars[pars$cond == "A", ]$kappa[i])
    B <- rmixture2p(N_obs, mu = 0, pMem = pars[pars$cond == "B", ]$pmem[i],
                    kappa = pars[pars$cond == "B", ]$kappa[i])
    data[[i]] <- c(A, B)
  }
  
  data <- data.frame(id = rep(1:N_subj, each = 2 * N_obs), 
                     y = unlist(data),
                     target = rep(0, 2 * N_subj * N_obs),
                     cond = rep(c("A","B"), each = N_obs))
  
  par_ml <- fit_mixtur(data, model = "2_component", unit = "radians",
                       response_var = "y", condition_var = "cond")
  names(par_ml) <- c('id','kappa','pmem','pguess', 'cond')
  par_ml <- par_ml[,c('id','cond', 'pmem','kappa')]
  
  fit_bmm <- list()
  fit_bmm[[1]] <- bmm(bmf(thetat ~ 1 + (cond || id),
                         kappa ~ 1 + (cond || id)),
                     data = data,
                     model = mixture2p('y'),
                     cores = cores,
                     backend = 'cmdstanr',
                     save_pars = save_pars(all = TRUE),
                     sample_prior = TRUE,
                     refresh = 0)
  
  fit_bmm[[2]] <- bmm(bmf(thetat ~ cond + (cond || id),
                          kappa ~ 1 + (cond || id)),
                      data = data,
                      model = mixture2p('y'),
                      prior = set_prior("normal(0, 1)", class = "b", nlpar = "thetat"),
                      cores = cores,
                      backend = 'cmdstanr',
                      save_pars = save_pars(all = TRUE),
                      sample_prior = TRUE,
                      refresh = 0)
      
  fit_bmm[[3]] <- bmm(bmf(thetat ~ 1 + (cond || id),
                          kappa ~ cond + (cond || id)),
                      data = data,
                      model = mixture2p('y'),
                      cores = cores,
                      backend = 'cmdstanr',
                      save_pars = save_pars(all = TRUE),
                      sample_prior = TRUE,
                      refresh = 0)
      
  fit_bmm[[4]] <- bmm(bmf(thetat ~ cond + (cond || id),
                      kappa ~ cond + (cond || id)),
                      data = data,
                      model = mixture2p('y'),
                      prior = set_prior("normal(0, 1)", class = "b", nlpar = "thetat"),
                      cores = cores,
                      backend = 'cmdstanr',
                      save_pars = save_pars(all = TRUE),
                      sample_prior = TRUE,
                      refresh = 0)
          
  par_bmm <- list()
  for (i in 1:length(fit_bmm)) {
    par_names <- grep('thetat_|kappa_', rownames(fixef(fit_bmm[[i]])), value = TRUE)
    all_pars <- c('kappa_Intercept','kappa_condB','thetat_Intercept','thetat_condB')
    fe_bmm <- fixef(fit_bmm[[i]])[par_names,1]
    missing <- setdiff(all_pars, par_names)
    fe_bmm[missing] <- 0
    fe_bmm <- fe_bmm[c('kappa_Intercept','kappa_condB','thetat_Intercept','thetat_condB')]
    
    
    re_bmm <- ranef(fit_bmm[[i]])$id[,1,]
    par_bmm_tmp <- as.data.frame(t(t(re_bmm) + fe_bmm))
    par_bmm_tmp$id <- 1:N_subj
    par_bmm_tmp$kappa_condB <- par_bmm_tmp$kappa_condB + par_bmm_tmp$kappa_Intercept
    par_bmm_tmp$thetat_condB <- par_bmm_tmp$thetat_condB + par_bmm_tmp$thetat_Intercept
    
    par_bmm_tmp <- data.frame(id = rep(1:N_subj, each = 2),
                              pmem = c(par_bmm_tmp$thetat_Intercept, par_bmm_tmp$thetat_condB),
                              kappa = c(par_bmm_tmp$kappa_Intercept, par_bmm_tmp$kappa_condB),
                              cond = rep(c("A","B"), each = N_subj))
    
    par_bmm_tmp$kappa <- exp(par_bmm_tmp$kappa)
    par_bmm_tmp$pmem <- exp(par_bmm_tmp$pmem)/(1+exp(par_bmm_tmp$pmem))
    attr(par_bmm_tmp, "formula") <- fit_bmm[[i]]$bmm$user_formula
    class(par_bmm_tmp) <- c('bmmpars', 'data.frame')
    
    par_bmm[[i]] <- par_bmm_tmp
  }
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


par_rec_stats_cond <- function(fits) {
  n_bmm <- length(fits$par_bmm)
  with(fits, {
    data.frame(
      formula = rep(c('',sapply(par_bmm, function(x) {
        paste0(capture.output(print(attr(x, 'formula'))), collapse=', ')
        })), each = 2),
      engine = c('ML','ML', rep('BMM', 2*n_bmm)),
      param = rep(c('pmem','kappa'), n_bmm + 1),
      r = c(cor.test(pars$pmem, par_ml$pmem)$estimate,
            cor.test(pars$kappa, par_ml$kappa)$estimate,
            sapply(par_bmm, function(x) {
              c(cor.test(pars$pmem, x$pmem)$estimate,
                cor.test(pars$kappa, x$kappa)$estimate)
              })),
      rmse = c(sqrt(mean((pars$pmem - par_ml$pmem)^2)),
               sqrt(mean((pars$kappa - par_ml$kappa)^2)),
               sapply(par_bmm, function(x) {
                 c(sqrt(mean((pars$pmem - x$pmem)^2)),
                   sqrt(mean((pars$kappa - x$kappa)^2)))
                 })),
      pop_par_diff = c(mean(pars$pmem) - mean(par_ml$pmem),
                       mean(pars$kappa) - mean(par_ml$kappa),
                       sapply(par_bmm, function(x) {
                         c(mean(pars$pmem) - mean(x$pmem),
                           mean(pars$kappa) - mean(x$kappa))
                         }))
    )
  })
}


print.bmmpars <- function(x) {
  cat('Formula: \n')
  print(attr(x, 'formula'))
  cat('\n\n')
  class(x) <- 'data.frame'
  print(x)
}
