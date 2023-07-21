#' This is the tutorial script for setting up the Bays et al. (2009) mixture model
#' for visual working memory tasks that use continuous report recall procedures.
#' 
#' In this script, you will see:
#'  1) how the model is set up using the brms package, 
#'  2) how a simple version of the model is estimates, and 
#'  3) how the model can be evaluated and results extracted and plotted.

# 0) R Setup: Packages & Data --------------------------------------------------
# start fresh
rm(list = ls()) # clean up work space
graphics.off()  # switch off graphics device

# load required packages
library(brms)       # for estimating the mixture model
library(tidyverse)
library(here)
library(bmm)

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

# load data file
data_Bays2009 <- read.table(here("data/Bays2009.txt"), header = T) %>% 
  dplyr::mutate(
    RespErr = wrap(RespErr),
    # create index variables for whether the lure exist at that setsize
    LureIdx1 = case_when(setsize >= 2 ~ 1, TRUE ~ 0),
    LureIdx2 = case_when(setsize >= 3 ~ 1, TRUE ~ 0),
    LureIdx3 = case_when(setsize >= 4 ~ 1, TRUE ~ 0),
    LureIdx4 = case_when(setsize >= 5 ~ 1, TRUE ~ 0),
    LureIdx5 = case_when(setsize >= 6 ~ 1, TRUE ~ 0),
    # recode lures to be relative to target
    Pos_Lure1 = wrap(RespErr-Pos_Lure1),
    Pos_Lure2 = wrap(RespErr-Pos_Lure2),
    Pos_Lure3 = wrap(RespErr-Pos_Lure3),
    Pos_Lure4 = wrap(RespErr-Pos_Lure4),
    Pos_Lure5 = wrap(RespErr-Pos_Lure5),
    # variable to include in formula as a correction to theta due to setsize
    inv_ss = 1/(setsize-1),
    inv_ss = ifelse(is.infinite(inv_ss), 1, inv_ss),
    setsize = as.factor(setsize)) %>% 
    select(subID,trial,setsize,RespErr,
         Pos_Lure1, Pos_Lure2, Pos_Lure3, Pos_Lure4, Pos_Lure5,
         LureIdx1, LureIdx2, LureIdx3, LureIdx4, LureIdx5, inv_ss)

data_Bays2009[is.na(data_Bays2009)] <- 0
head(data_Bays2009[sample(1:nrow(data_Bays2009),6),])

###############################################################################!
# 2) BRMS fit ------------------------------------------------------------------
###############################################################################!


# create mixture of von Mises distributions
Bays_mixModel <- mixture(von_mises(link = "identity"),
                         von_mises(link = "identity"),
                         von_mises(link = "identity"),
                         von_mises(link = "identity"),
                         von_mises(link = "identity"),
                         von_mises(link = "identity"),
                         von_mises(link = "identity"),
                         order = "none")

# set up mixture model
Bays_mixModel_formula <- bf(RespErr ~ 1,
                            # fix kappa over memory distributions
                            nlf(kappa1 ~ kappa),     # target distribution
                            nlf(kappa2 ~ kappa),     # non-target
                            nlf(kappa3 ~ kappa),     # non-target
                            nlf(kappa4 ~ kappa),     # non-target
                            nlf(kappa5 ~ kappa),     # non-target
                            nlf(kappa6 ~ kappa),     # non-target
                            # kappa for guessing distribution will be fixed using priors
                            kappa7 ~ 1,             # uniform  
                            # specify mixing distributions for distinct item categories
                            nlf(theta1 ~ thetat),   # p_mem
                            nlf(theta2 ~ LureIdx1*(thetant + log(inv_ss)) + (1-LureIdx1)*(-100)),  # p_intrusion
                            nlf(theta3 ~ LureIdx2*(thetant + log(inv_ss)) + (1-LureIdx2)*(-100)),  # p_intrusion
                            nlf(theta4 ~ LureIdx3*(thetant + log(inv_ss)) + (1-LureIdx3)*(-100)),  # p_intrusion
                            nlf(theta5 ~ LureIdx4*(thetant + log(inv_ss)) + (1-LureIdx4)*(-100)),  # p_intrusion
                            nlf(theta6 ~ LureIdx5*(thetant + log(inv_ss)) + (1-LureIdx5)*(-100)),  # p_intrusion
                            # target & guessing distribution will be centered using priors
                            mu1 ~ 1, # fixed intercept constrained using priors
                            mu7 ~ 1, # fixed intercept constrained using priors
                            # center non-target distribution on data specified locations
                            nlf(mu2 ~ Pos_Lure1),           # center non-target
                            nlf(mu3 ~ Pos_Lure2),           # center non-target
                            nlf(mu4 ~ Pos_Lure3),           # center non-target
                            nlf(mu5 ~ Pos_Lure4),           # center non-target
                            nlf(mu6 ~ Pos_Lure5),           # center non-target
                            # now predict parameters of interest
                            kappa ~ 0 + setsize,     # fixed intercept for precision of memory distributions
                            thetat ~ 0 + setsize,    # fixed intercept for p_mem
                            thetant ~ 0 + setsize,   # fixed intercept for p_intrusion
                            # thetag ~ 1,    # fixed intercept for p_guess
                            # for brms to process this formula correclty, set non-linear to TRUE
                            nl = TRUE)


# check default priors
get_prior(Bays_mixModel_formula, data_Bays2009, Bays_mixModel)

# constrain priors to identify the model
Bays_mixModel_priors <- 
  # first we center the target and guessing distribution to zero
  prior(constant(0), class = Intercept, dpar = "mu1") + 
  prior(constant(0), class = Intercept, dpar = "mu7") +
  # next, we set the guessing distribution to be uniform, kappa -> 0
  prior(constant(-100), class = Intercept, dpar = "kappa7") +
  # next, we set reasonable priors for the to be estimated distributions
  prior(normal(5.0, 0.8), class = b, nlpar = "kappa") +
  prior(logistic(0, 1), class = b, nlpar = "thetat") +
  prior(logistic(0, 1), class = b, nlpar = "thetant") +
  prior(constant(-100), class = b, coef="setsize1", nlpar="thetant")

# if the model has been already estimated, load the results, otherwise estimate it
filename <- 'output/fit_bays2009_3p_model.RData'
if (!file.exists(here(filename))) {
  # fit the mixture model using brms
  fit_Bays_mixMod <- brm(
    # include model information
    formula = Bays_mixModel_formula, # specify formula for mixture model
    data    = data_Bays2009, # specify data used to estimate the mixture model
    family  = Bays_mixModel, # call the defined mixture family
    prior   = Bays_mixModel_priors, # use the used defined priors,
    
    # save settings
    sample_prior = TRUE,
    save_pars = save_pars(all = TRUE),
    
    # add brms settings
    warmup = warmup_samples,
    iter = warmup_samples + postwarmup_samples, 
    chains = nChains,
    
    # control commands for the sampler
    control = list(adapt_delta = adapt_delta, 
                   max_treedepth = max_treedepth)
  )
  
  # save results into file
  save(fit_Bays_mixMod, 
       file = here(filename),
       compress = "xz")
} else {
  # load results file
  load(file = here(filename))
}

#############################################################################!
# 2) Model evaluation                                                    ####
#############################################################################!

## 2.1) fit & summary ----------------------------------------------------------
# plot the posterior predictive check to evaluate overall model fit
pp_check(fit_Bays_mixMod)

# print results summary. There is 1 divergent transition, but we will ignore it
# for the purposes of this illustration. For real analyses, follow the suggestions
# to remove this issue
summary(fit_Bays_mixMod)

## 2.2) Extract parameter estimates --------------------------------------------
# extract the fixed effects from the model and determine the rows that contain
# the relevant parameter estimates
fixedEff <- fixef(fit_Bays_mixMod)
thetat <- fixedEff[grepl("thetat",rownames(fixedEff)),]
thetant <- fixedEff[grepl("thetant",rownames(fixedEff)),]
kappa <- fixedEff[grepl("kappa_",rownames(fixedEff)),]

# transform parameters because brms uses special link functions
kappa <- exp(kappa)
sd <- k2sd(kappa[,1]) 
pmem <- exp(thetat)/(exp(thetat)+exp(thetant)+exp(0))
pnt <- exp(thetant)/(exp(thetat)+exp(thetant)+exp(0))
pg <- exp(0)/(exp(thetat)+exp(thetant)+exp(0))

# print parameter estimates
kappa
sd #in radians
round(pmem,3)
round(pnt,3)
round(pg, 3)



## 2.3) plot parameter estimates -----------------------------------------------
# define defaults for clean ggplots
clean_plot <- theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    axis.line.x = element_line(color = 'black'),
                    axis.line.y = element_line(color = 'black'),
                    legend.key = element_rect(fill = 'white'),
                    text = element_text(size = 15),
                    line = element_line(linewidth = 1),
                    axis.ticks = element_line(linewidth = 1))

# extract brms draws from the posterior distribution
fixedFX_draws <- fit_Bays_mixMod %>% 
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
         pmem = exp(thetat)/(exp(thetat)+exp(thetant)+exp(0)),
         pnt = exp(thetant)/(exp(thetat)+exp(thetant)+exp(0)),
         pg = 1-pmem-pnt,
         setsize = str_remove_all(cond,"setsize"))

# results from the published paper
results_Bays2009 <- data.frame(
  setsize = as.character(c(1,2,4,6)),
  pnt = c(0,0.03,0.115,0.28),
  pg = c(0.01,0.05,0.16,0.14),
  sd = c(13.5,18,22.5,24.5)
)
results_Bays2009$pmem <- 1-results_Bays2009$pnt-results_Bays2009$pg

# plot sd results
(sd_plot <- ggplot(fixedFX_draws, aes(x = setsize, y = sd)) +
    geom_half_violin(position = position_nudge(x = .10, y = 0), side = "r", fill = "darkgrey", color = NA,
                     alpha = 0.9, scale = "width") +
    stat_summary(geom = "pointrange", fun.data = mode_hdci, color = "black",
                 size = 0.7, linewidth = 0.8,
                 position = position_dodge(0.25)) +
    geom_point(data = results_Bays2009,
               aes(x = setsize, y = sd), color = "black",
               shape = "diamond", size = 4,
               position = position_nudge(x = .25, y = 0)) +
    scale_fill_grey(start = 0, end = .8) +
    scale_color_grey(start = 0, end = .8) +
    labs(x = "Set Size", y = "Memory imprecision (SD)", title = "A") +
    clean_plot)

# plot pnt results
(pnt_plot <- ggplot(fixedFX_draws, aes(x = setsize, y = pnt)) +
    geom_half_violin(position = position_nudge(x = .10, y = 0), side = "r", fill = "darkgrey", color = NA,
                     alpha = 0.9, scale = "width") +
    stat_summary(geom = "pointrange", fun.data = mean_hdci, color = "black",
                 size = 0.7, linewidth = 0.8,
                 position = position_dodge(0.25)) +
    geom_point(data = results_Bays2009,
               aes(x = setsize, y = pnt), color = "black",
               shape = "diamond", size = 4,
               position = position_nudge(x = .25, y = 0)) +
    scale_fill_grey(start = 0, end = .8) +
    scale_color_grey(start = 0, end = .8) +
    labs(x = "Set Size", y = "Non-target responses", title = "B") +
    clean_plot)

# plot pg results
(pg_plot <- ggplot(fixedFX_draws, aes(x = setsize, y = pg)) +
    geom_half_violin(position = position_nudge(x = .10, y = 0), side = "r", fill = "darkgrey", color = NA,
                     alpha = 0.9, scale = "width") +
    stat_summary(geom = "pointrange", fun.data = mode_hdci, color = "black",
                 size = 0.7, linewidth = 0.8,
                 position = position_dodge(0.25)) +
    geom_point(data = results_Bays2009,
               aes(x = setsize, y = pg), color = "black",
               shape = "diamond", size = 4,
               position = position_nudge(x = .25, y = 0)) +
    scale_fill_grey(start = 0, end = .8) +
    scale_color_grey(start = 0, end = .8) +
    labs(x = "Set Size", y = "Random responses", title = "C") +
    clean_plot)

(sd_plot | pnt_plot | pg_plot)

ggsave(here("figures/plot_Bays2009.jpeg"), width = 9, height = 3)
