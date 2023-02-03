#' This is the code for the unconstrained version of Example 4 from the paper We
#' fit the 2 parameter model without constraining the parameters to be monotonic
#'
#' 

#############################################################################!
# 0) R Setup                                                             ####
#############################################################################!

pacman::p_load(here, brms, tidyverse, tidybayes, patchwork, gghalves)
pacman::p_load_gh("venpopov/bmm")
source(here('functions/clean_plot.R'))

data_popov <- read.csv(here('data/popov_unpublished.csv'))
data_popov$setsize <- as.factor(data_popov$setsize)

#############################################################################!
# 1) Model specification and estimation                                  ####
#############################################################################!

# code contrasts such that the effect of each size is relative to the effect of 
# the setsize before
contr_ordin <- matrix(0, nrow=8, ncol=7)
contr_ordin[lower.tri(contr_ordin)] <- -1
colnames(contr_ordin) <- 2:8
contrasts(data_popov$setsize) <- contr_ordin
contrasts(data_popov$setsize)

ff <- bf(anglediff ~ 1,
         kappa ~ setsize,
         thetat ~ setsize)

pr <- prior_('normal(0.0,0.8)', class='b', nlpar='kappa', lb=0) +
  prior_('logistic(0,1)', class='b', nlpar='thetat',lb=0)



filename <- 'output/fit_popov_unpublished_monotonic_constraints_2p.RData'
if (!file.exists(here(filename))) {
  
  fit_popov <- fit_model(formula = ff,
                         data = data_popov,
                         model_type = "2p",
                         parallel = TRUE,
                         warmup = 500,
                         iter = 1000,
                         prior=pr)
  
  save(list=ls(), file=here(filename))  
} else load(here(filename))

#############################################################################!
# 2) Model evaluation                                                    ####
#############################################################################!


summary(fit_popov)

fixedEff <- fixef(fit_popov)
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
fixedFX_draws <- fit_popov %>% 
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
    coord_cartesian(ylim=c(0,0.8)) +
    labs(x = "Set Size", y = "Random responses", title = "D") +
    clean_plot())

(sd_plot | pg_plot)

ggsave(here("figures/plot_popov_unpublished_monotonic_constraints.jpeg"), width = 6, height = 3)
