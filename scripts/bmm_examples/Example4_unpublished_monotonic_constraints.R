#' This is the code for the unconstrained version of Example 4 from the paper We
#' fit the 2 parameter model constraining the parameters to be monotonic

#############################################################################!
# 0) R Setup                                                             ####
#############################################################################!

pacman::p_load(here, brms, tidyverse, tidybayes, patchwork, gghalves)
pacman::p_load_gh("venpopov/bmm")

source(here('functions/clean_plot.R'))

# load missing output files
source(here("scripts","LoadResultsFiles.R"))

# specify the number of samples to run for warm up & after warm up
warmup_samples <- 1000
postwarmup_samples <- 1000

# specify the number of chains
nChains <- 4

#' if the number of user defined chains is larger than the number of cores 
#' on the system than estimate as many chains as there are cores on the system
if (nChains >  parallel::detectCores()) {
  nChains <-  parallel::detectCores()
}

# set brms controls to prevent divergent transitions & improve convergence
adapt_delta <- .95
max_treedepth <- 10

data_popov <- read.csv(here('data/popov_unpublished.csv'))
data_popov$setsize <- as.factor(data_popov$setsize)

#############################################################################!
# 1) Model specification and estimation                                  ####
#############################################################################!

model_2p <- mixture2p(resp_error = "anglediff")

# code contrasts such that the effect of each size is relative to the effect of 
# the setsize before
contr_ordin <- matrix(0, nrow = 8, ncol = 7)
contr_ordin[lower.tri(contr_ordin)] <- -1
colnames(contr_ordin) <- 2:8
contrasts(data_popov$setsize) <- contr_ordin
contrasts(data_popov$setsize)

ff <- bmf(kappa ~ 1 + setsize + (1 + setsize |subject),
          thetat ~ 1 + setsize + (1 + setsize |subject))

pr <- prior_('normal(0.0,0.8)', class = 'b', nlpar = 'kappa', lb = 0) +
  prior_('logistic(0,1)', class = 'b', nlpar = 'thetat',lb = 0)

filename <- here('output','fit_E4_unpublished_monotonic_constraints_2p')

fit_popov_prior <- bmm(model = model_2p,
                       formula = ff,
                       data = data_popov,
                       prior = pr,
                       
                       # save settings
                       sample_prior = TRUE,
                       save_pars = save_pars(all = TRUE),
                       
                       # add brms settings
                       warmup = warmup_samples,
                       iter = warmup_samples + postwarmup_samples, 
                       chains = nChains,
                       cores = parallel::detectCores(),
                       
                       # control commands for the sampler
                       control = list(adapt_delta = adapt_delta, 
                                      max_treedepth = max_treedepth),
                       
                       # save results to file
                       file = filename)

fit_popov <- readRDS(here("output","fit_popov_unpublished_noconstraints_2p.rds"))

bf <- bayes_factor(fit_popov,fit_popov_prior, repetitions = 10, cores = 4)
hist(log(bf$bf))

#############################################################################!
# 2) Model evaluation                                                    ####
#############################################################################!

pp_check(fit_popov_prior, group = "setsize", type = "dens_overlay_grouped")

summary(fit_popov_prior)

fixedEff <- fixef(fit_popov_prior)
kappa <- fixedEff[grepl('kappa_', row.names(fixedEff)),]
kappa[2:nrow(kappa),] <- -kappa[2:nrow(kappa),]
kappa <- apply(kappa,2,cumsum)
sd <- apply(kappa, 2, function(x) k2sd(exp(x)))*180/pi
thetat <- fixedEff[grepl('thetat_', row.names(fixedEff)),]
thetat[2:nrow(thetat),] <- -thetat[2:nrow(thetat),]
thetat <- apply(thetat,2,cumsum)
pmem <- exp(thetat)/(exp(thetat)+exp(0))

pars <- data.frame(par=rep(c('sd','pmem'), each=8),
                   setsize=rep(1:8,2),
                   value= c(sd[,1],pmem[,1]),
                   Q2.5 =c(sd[,3],pmem[,3]),
                   Q97.5 = c(sd[,4],pmem[,4]))

ggplot(pars, aes(setsize, value)) +
  geom_point() +
  geom_line() +
  facet_wrap(~par, scales='free')


# extract brms draws from the posterior distribution
fixedFX_draws <- fit_popov_prior %>% 
  tidy_draws() %>%
  select(starts_with("b_"),.chain,.iteration,.draw) %>%  
  pivot_longer(cols = starts_with("b_"),
               names_to = "modelPar",
               values_to = "postSample") %>% 
  mutate(par = str_split_i(modelPar,"_",2),
         cond = str_split_i(modelPar,"_",3)) %>% 
  select(-modelPar) %>% 
  filter(par %in% c('kappa','thetat')) %>%
  mutate(postSample = ifelse(cond=="Intercept",postSample,-postSample)) %>% 
  group_by(.chain,.iteration,.draw,par) %>% 
  mutate(postSample = cumsum(postSample)) %>% 
  spread(par, postSample) %>% 
  mutate(sd = k2sd(exp(kappa))/pi * 180,
         pmem = exp(thetat)/(exp(thetat)+exp(0)),
         pg = 1-pmem,
         setsize = str_remove_all(cond,"setsize"),
         setsize = ifelse(setsize=="Intercept",1,setsize))

# plot sd results
(sd_plot <- ggplot(fixedFX_draws, aes(x = setsize, y = sd)) +
    geom_half_violin(position = position_nudge(x = .05, y = 0), side = "r", fill = "darkgrey", color = NA,
                     alpha = 0.9, scale = "width") +
    stat_summary(geom = "pointrange", fun.data = mode_hdi, color = "black",
                 size = 0.3, linewidth = 0.8,
                 position = position_dodge(0.1)) +
    scale_fill_grey(start = 0, end = .8) +
    scale_color_grey(start = 0, end = .8) +
    coord_cartesian(ylim=c(10,45)) +
    labs(x = "Set Size", y = "Memory imprecision (SD)", title = "C") +
    clean_plot())

# plot pnt results
(pg_plot <- ggplot(fixedFX_draws, aes(x = setsize, y = pg)) +
    geom_half_violin(position = position_nudge(x = .05, y = 0), side = "r", fill = "darkgrey", color = NA,
                     alpha = 0.9, scale = "width") +
    stat_summary(geom = "pointrange", fun.data = mode_hdi, color = "black",
                 size = 0.3, linewidth = 0.8,
                 position = position_dodge(0.1)) +
    scale_fill_grey(start = 0, end = .8) +
    scale_color_grey(start = 0, end = .8) +
    coord_cartesian(ylim = c(0,0.8)) +
    labs(x = "Set Size", y = "Random responses", title = "D") +
    clean_plot())

(sd_plot | pg_plot)

ggsave(here("figures/plot_popov_unpublished_monotonic_constraints.jpeg"), width = 6, height = 3)
