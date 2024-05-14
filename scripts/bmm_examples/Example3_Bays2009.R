#' This is a tutorial for estimating the Bays et al (2009) three-parameter
#' mixture model for visual working memory tasks that use continuous report
#' recall procedures. We will use the `bmm` package, which provides wrapper
#' functions around brms. The `bmm` package makes it easy to specify the three
#' parameter model.
#' 
#' 
# start fresh
rm(list = ls()) # clean up work space
graphics.off()  # switch off graphics device

pacman::p_load(here, brms, tidyverse, tidybayes, patchwork, gghalves)
pacman::p_load_gh("venpopov/bmm")

# load function to clean up plots
source(here("functions","clean_plot.R"))

# load missing output files
source(here("scripts","LoadResultsFiles.R"))

# Set up parallel sampling of mcmc chains
options(mc.cores =  parallel::detectCores())

# specify the number of samples to run for warm up & after warm up
warmup_samples <- 2000
postwarmup_samples <- 2000

# specify the number of chains
nChains <- 6

#' if the number of user defined chains is larger than the number of cores 
#' on the system than estimate as many chains as there are cores on the system
if (nChains >  parallel::detectCores()) {
  nChains <-  parallel::detectCores()
}

# set brms controls to prevent divergent transitions & improve convergence
adapt_delta <- .95
max_treedepth <- 10

#############################################################################!
# 0) R Setup                                                             ####
#############################################################################!

data_bays_2009 <- read.table(here("data/Bays2009.txt"), header = T)

# in the data, Pos_Lure1 to Pos_Lure5 contain the error relative to each lure
# position. To estimate the model, we need the actual position of the lure 
# relative to the target. Since RespErr is error relative to the target, we need
# to subtract them. We also use the wrap function from the bmm package to ensure
# that the results are between -pi and pi on the circle
data_bays_2009 <- data_bays_2009 %>% 
  mutate(RespErr = wrap(RespErr),
         col_lure1 = wrap(RespErr - Pos_Lure1),
         col_lure2 = wrap(RespErr - Pos_Lure2),
         col_lure3 = wrap(RespErr - Pos_Lure3),
         col_lure4 = wrap(RespErr - Pos_Lure4),
         col_lure5 = wrap(RespErr - Pos_Lure5),
         setsize = as.factor(setsize)) %>% 
  select(subID, setsize, RespErr, col_lure1:col_lure5) %>% 
  group_by(subID,setsize) %>% 
  mutate(trialID = 1:n()) %>%  ungroup()

data_bays_2009[is.na(data_bays_2009)] <- 0
head(data_bays_2009 %>% arrange(subID, trialID))

#############################################################################!
# 1) Model specification and estimation                                  ####
#############################################################################!

# model object
model_bays_2009 <- mixture3p(resp_error = "RespErr",
                             nt_features = "col_lure",
                             set_size = "setsize",
                             regex = TRUE)

# formula
formula_bays_2009 <- bmf(kappa ~ 0 + setsize + (0 + setsize || subID),
                         thetat ~ 0 + setsize + (0 + setsize || subID),
                         thetant ~ 0 + setsize + (0 + setsize || subID))

# if the model has been already estimated, load the results, otherwise estimate it
filename <- here('output','fit_bays2009_3p_model')

fit_bays_2009 <- bmm(model = model_bays_2009,
                     formula = formula_bays_2009, 
                     data = data_bays_2009, 
                     
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

#############################################################################!
# 2) Model evaluation                                                    ####
#############################################################################!

## 2.1) fit & summary ----------------------------------------------------------
# plot the posterior predictive check to evaluate overall model fit
pp_check(fit_bays_2009)

# print resultformula# print results summary. There is 1 divergent transition, but we will ignore it
# for the purposes of this illustration. For real analyses, follow the suggestions
# to remove this issue
summary(fit_bays_2009)

## 2.2) Extract parameter estimates --------------------------------------------
# extract the fixed effects from the model and determine the rows that contain
# the relevant parameter estimates
fixedEff <- fixef(fit_bays_2009)
thetat <- fixedEff[grepl("thetat",rownames(fixedEff)),]
thetant <- fixedEff[grepl("thetant",rownames(fixedEff)),]
kappa <- fixedEff[grepl("kappa_",rownames(fixedEff)),]

# transform parameters because brms uses special link functions
kappa <- exp(kappa)
sd <- k2sd(kappa[,1]) 
pmem <- exp(thetat) / (exp(thetat) + exp(thetant) + exp(0))
pnt <- exp(thetant) / (exp(thetat) + exp(thetant) + exp(0))
pg <- exp(0) / (exp(thetat) + exp(thetant) + exp(0))

# print parameter estimates
kappa
sd #in radians
round(pmem,3)
round(pnt,3)
round(pg, 3)



## 2.3) plot parameter estimates -----------------------------------------------
# extract brms draws from the posterior distribution
fixedFX_draws <- fit_bays_2009 %>% 
  tidy_draws() %>%
  select(starts_with("b_"),.chain,.iteration,.draw) %>% 
  pivot_longer(cols = starts_with("b_"),
               names_to = "modelPar",
               values_to = "postSample") %>% 
  mutate(par = str_split_i(modelPar,"_",2),
         cond = str_split_i(modelPar,"_",3)) %>% 
  select(-modelPar) %>% 
  filter(par %in% c('kappa','thetat','thetant')) %>%
  spread(par, postSample) %>% 
  mutate(sd = k2sd(exp(kappa))/pi * 180,
         pmem = exp(thetat) / (exp(thetat) + exp(thetant) + exp(0)),
         pnt = exp(thetant) / (exp(thetat) + exp(thetant) + exp(0)),
         pg = 1 - pmem - pnt,
         setsize = str_remove_all(cond,"setsize"))

# results from the published paper
results_Bays2009 <- data.frame(
  setsize = as.character(c(1,2,4,6)),
  pnt = c(0,0.03,0.115,0.28),
  pg = c(0.01,0.05,0.16,0.14),
  sd = c(13.5,18,22.5,24.5)
)
results_Bays2009$pmem <- 1 - results_Bays2009$pnt - results_Bays2009$pg

# plot sd results
(sd_plot <- ggplot(fixedFX_draws, aes(x = setsize, y = sd)) +
    geom_half_violin(position = position_nudge(x = .05, y = 0), side = "r", fill = "darkgrey", color = NA,
                     alpha = 0.9, scale = "width") +
    stat_summary(geom = "pointrange", fun.data = mode_hdci, color = "black",
                 size = 0.3, linewidth = 0.8,
                 position = position_dodge(0.1)) +
    geom_point(data = results_Bays2009,
               aes(x = setsize, y = sd), color = "black",
               shape = "diamond", size = 2.5,
               position = position_nudge(x = .1, y = 0)) +
    scale_fill_grey(start = 0, end = .8) +
    scale_color_grey(start = 0, end = .8) +
    coord_cartesian(ylim = c(0,35)) +
    labs(x = "Set Size", y = "Memory imprecision (SD)", title = "A") +
    clean_plot())

# plot pnt results
(pnt_plot <- ggplot(fixedFX_draws, aes(x = setsize, y = pnt)) +
    geom_half_violin(position = position_nudge(x = .05, y = 0), side = "r", fill = "darkgrey", color = NA,
                     alpha = 0.9, scale = "width") +
    stat_summary(geom = "pointrange", fun.data = mode_hdci, color = "black",
                 size = 0.3, linewidth = 0.8,
                 position = position_dodge(0.1)) +
    geom_point(data = results_Bays2009,
               aes(x = setsize, y = pnt), color = "black",
               shape = "diamond", size = 2.5,
               position = position_nudge(x = .1, y = 0)) +
    scale_fill_grey(start = 0, end = .8) +
    scale_color_grey(start = 0, end = .8) +
    coord_cartesian(ylim = c(0,.5)) +
    labs(x = "Set Size", y = "Non-target responses", title = "B") +
    clean_plot())

# plot pg results
(pg_plot <- ggplot(fixedFX_draws, aes(x = setsize, y = pg)) +
    geom_half_violin(position = position_nudge(x = .05, y = 0), side = "r", fill = "darkgrey", color = NA,
                     alpha = 0.9, scale = "width") +
    stat_summary(geom = "pointrange", fun.data = mode_hdci, color = "black",
                 size = 0.3, linewidth = 0.8,
                 position = position_dodge(0.1)) +
    geom_point(data = results_Bays2009,
               aes(x = setsize, y = pg), color = "black",
               shape = "diamond", size = 2.5,
               position = position_nudge(x = .1, y = 0)) +
    scale_fill_grey(start = 0, end = .8) +
    scale_color_grey(start = 0, end = .8) +
    coord_cartesian(ylim = c(0,.5)) +
    labs(x = "Set Size", y = "Random responses", title = "C") +
    clean_plot())

(sd_plot | pnt_plot | pg_plot)

ggsave(here("figures/plot_Bays2009.jpeg"), width = 9, height = 3)
