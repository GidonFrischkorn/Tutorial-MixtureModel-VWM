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

ff <- bf(anglediff ~ 1,
         kappa ~ 0+setsize + (1|subject),
         thetat ~ 0+setsize + (1|subject))



filename <- 'output/fit_popov_unpublished_noconstraints_2p.RData'
if (!file.exists(here(filename))) {
  
  fit_popov <- fit_model(formula = ff,
                         data = data_popov,
                         model_type = "2p",
                         parallel = TRUE,
                         warmup = 1000,
                         iter = 2000)
  
  save(list=ls(), file=here(filename))  
} else load(here(filename))

#############################################################################!
# 2) Model evaluation                                                    ####
#############################################################################!


summary(fit_popov)

fixedEff <- fixef(fit_popov)
kappa <- fixedEff[grepl('kappa_', row.names(fixedEff)),]
sd <- apply(kappa, 2, k2sd)*180/pi
thetat <- fixedEff[grepl('thetat_', row.names(fixedEff)),]
pmem <- exp(thetat)/(exp(thetat)+exp(0))

pars <- data.frame(par=rep(c('sd','pmem'), each=8),
                   setsize=rep(1:8,2),
                   value= c(sd[,1],pmem[,1]),
                   Q2.5 =c(sd[,3],pmem[,3]),
                   Q97.5 = c(sd[,4],pmem[,4]))

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
  spread(par, postSample) %>% 
  mutate(sd = k2sd(exp(kappa))/pi * 180,
         pmem = exp(thetat)/(exp(thetat)+exp(0)),
         pg = 1-pmem,
         setsize = str_remove_all(cond,"setsize"))

# plot sd results
(sd_plot <- ggplot(fixedFX_draws, aes(x = setsize, y = sd)) +
    geom_half_violin(position = position_nudge(x = .05, y = 0), side = "r", fill = "darkgrey", color = NA,
                     alpha = 0.9, scale = "width") +
    stat_summary(geom = "pointrange", fun.data = mode_hdi, color = "black",
                 size = 0.3, linewidth = 0.8,
                 position = position_dodge(0.1)) +
    scale_fill_grey(start = 0, end = .8) +
    scale_color_grey(start = 0, end = .8) +
    labs(x = "Set Size", y = "Memory imprecision (SD)", title = "A") +
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
    labs(x = "Set Size", y = "Random responses", title = "B") +
    clean_plot())

(sd_plot | pg_plot)

ggsave(here("figures/plot_popov_unpublished_unconstrained.jpeg"), width = 6, height = 3)
