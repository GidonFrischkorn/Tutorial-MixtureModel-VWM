#' This is the tutorial script for setting up the Zhang & Luck (2008) mixture model
#' for visual working memory tasks that use continuous report recall procedures.
#' 
#' In this script, you will see:
#'  1) how the model is set up using the brms package, 
#'  2) how a simple version of the model is estimates, and 
#'  3) how the model can be evaluated and results extracted and plotted.

###############################################################################!
# 0) R Setup: Packages & Data --------------------------------------------------
###############################################################################!

# start fresh
rm(list = ls()) # clean up work space
graphics.off()  # switch off graphics device

# load required packages
pacman::p_load(here, brms, tidyverse, tidybayes, patchwork, gghalves, bmm)

# load function to clean up plots
source(here("functions","clean_plot.R"))

# load missing output files
source(here("scripts","LoadResultsFiles.R"))

# Set up parallel sampling of mcmc chains
options(mc.cores =  parallel::detectCores())

# specify the number of samples to run for warm up & after warm up
warmup_samples <- 3000
postwarmup_samples <- 3000

# specify the number of chains
nChains <- 6

#' if the number of user defined chains is larger than the number of cores 
#' on the system than estimate as many chains as there are cores on the system
if (nChains >  parallel::detectCores()) {
  nChains <-  parallel::detectCores()
}

# set brms controls to prevent divergent transitions & improve convergence
adapt_delta <- .99
max_treedepth <- 12

###############################################################################!
# 1) Read: Raw Data ------------------------------------------------------------
###############################################################################!

# load data file
data_ZL2008 <- read.table(here("data/Zhang&Luck2008.txt"), header = T) %>% 
  dplyr::mutate(setsize = as.factor(setsize),
                # wrap cases smaller than -pi,
                # or larger than pi around the circle
                RespErr = bmm::wrap(RespErr))

# have a look at the data and included variables
head(data_ZL2008)

###############################################################################!
# 2) Model Setup ---------------------------------------------------------------
###############################################################################!

#' To reproduce the results from Zhang & Luck (2008), we include setsize
#' as a within subject predictor.
#' Additionally we specify the formula so that we directly get the estimates
#' for each level of setsize: 0 + setsize
#' Finally, we implement the hierarchical structure by allowing that the setsize
#' effect varies over each subject: (0 + setsize || subID)
#' This is done for both kappa, the precision of the memory distribution,
#' and thetat, the mixing distribution for target responses, essentially estimating pMem.
ZL_mixFormula_bmm <- bf(
  # Initializing the dependent variable
  RespErr ~ 1,   
  # estimating fixed intercept & random intercept for kappa of the first von Mises
  kappa ~ 0 + setsize + (0 + setsize || subID), 
  # estimating fixed intercept & random intercept for the mixing proportion 
  # for the memory/target distribution (i.e., p_mem)
  thetat ~ 0 + setsize + (0 + setsize || subID)
)


###############################################################################!
# 3) Model estimation ----------------------------------------------------------
###############################################################################!

# fit mixture model if there is not already a results file stored
if (!file.exists(here("output/fit_E1_ZL2008.RData"))) {
  # using the fit_model function, we pass the bmm formula for the mixture model
  # prior specification and other constraints are taken care of within this function
  # additionally any further arguments for brms can be passed as well
  fit_ZL_mixModel <- fit_model(
    formula = ZL_mixFormula_bmm, # specify formula for mixture model
    data    = data_ZL2008,   # specify data used to estimate the mixture model
    model_type = "2p", # select the two-parameter model for fitting
    
    # save all potentially relevant information
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
  save(fit_ZL_mixModel, 
       file = here("output/fit_E1_ZL2008.RData"),
       compress = "xz")
  
} else {
  # load results file
  load(file = here("output/fit_E1_ZL2008.RData"))
}

###############################################################################!
# 4) Model evaluation ----------------------------------------------------------
###############################################################################!

## 4.1) fit & summary ----------------------------------------------------------
# plot the posterior predictive check to evaluate overall model fit
pp_check(fit_ZL_mixModel)

# print results summary
summary(fit_ZL_mixModel)

## 4.2) extract parameter estimates --------------------------------------------

# extract the fixed effects from the model
fixedEff <- fixef(fit_ZL_mixModel)

# determine the rows that contain the relevant parameter estimates
theta_cols <- grepl("theta",rownames(fixedEff))
kappa_cols <- grepl("kappa1",rownames(fixedEff))

# extract kappa estimates
kappa_fixedFX <- fixedEff[kappa_cols,]

# convert kappa estimates to absolute scale (radians)
kappa_fixedFX <- exp(kappa_fixedFX)

# extract theta estimates
theta_fixedFX <- fixedEff[theta_cols,]

# convert theta estimates into pMem estimates
p_Mem_fixedFX <- gtools::inv.logit(theta_fixedFX)

# print out parameter estimates
kappa_fixedFX
p_Mem_fixedFX

## 4.3) plot parameter estimates -----------------------------------------------
fixedFX_draws <- fit_ZL_mixModel %>% 
  tidy_draws() %>%
  select(starts_with("b_"),.chain,.iteration,.draw) %>% 
  pivot_longer(cols = starts_with("b_"),
               names_to = "modelPar",
               values_to = "postSample") %>% 
  mutate(par = str_split_i(modelPar,"_",2),
         cond = str_split_i(modelPar,"_",3)) %>% 
  select(-modelPar) %>% 
  filter(par == "kappa1" | par == "theta1") %>% 
  mutate(postSample_abs = case_when(par == "kappa1" ~ (sqrt(1/exp(postSample))/pi) * 180,
                                    par == "theta1" ~ inv_logit_scaled(postSample)),
         cond = str_remove_all(cond,"setsize"))

results_ZL2008 <- data.frame(
  cond = as.character(c(1,2,3,6)),
  pMem = c(.99,.95,.83,.38),
  sdMen = c(13.9,19.4,21.9,22.3)
)

# plot kappa results
kappa_plot <- ggplot(data = fixedFX_draws %>% filter(par == "kappa1"),
                     aes(x = cond, y = postSample_abs)) +
  coord_cartesian(ylim = c(5,45)) +
  geom_half_violin(position = position_nudge(x = .05, y = 0), side = "r", fill = "darkgrey", color = NA,
                   alpha = 0.9, scale = "width") +
  stat_summary(geom = "pointrange", fun.data = mean_hdci, color = "black",
               size = 0.7, linewidth = 0.8,
               position = position_dodge(0.2)) +
  geom_point(data = results_ZL2008,
             aes(y = sdMen, x = cond), color = "black",
             shape = "diamond", size = 4,
             position = position_nudge(x = .15, y = 0)) +
  scale_fill_grey(start = 0, end = .8) +
  scale_color_grey(start = 0, end = .8) +
  labs(x = "Set Size", y = "Memory imprecision (SD)", title = "B") +
  clean_plot()
kappa_plot

# plot pMem results
pMem_plot <- ggplot(data = fixedFX_draws %>% filter(par == "theta1"),
                    aes(x = cond, y = postSample_abs)) +
  coord_cartesian(ylim = c(0.2,1.05)) +
  geom_half_violin(position = position_nudge(x = .05, y = 0), side = "r", fill = "darkgrey", color = NA,
                   alpha = 0.9, scale = "width") +
  stat_summary(geom = "pointrange", fun.data = mean_hdci, color = "black",
               size = 0.7, linewidth = 0.8,
               position = position_dodge(0.2)) +
  geom_point(data = results_ZL2008,
             aes(y = pMem, x = cond), color = "black",
             shape = "diamond", size = 4,
             position = position_nudge(x = .15, y = 0)) +
  scale_fill_grey(start = 0, end = .8) +
  scale_color_grey(start = 0, end = .8) +
  labs(x = "Set Size", y = expression(P[mem]), title = "B") +
  clean_plot()
pMem_plot

# patch plots together
joint_plot <- (pMem_plot | kappa_plot)

# show joint plot
joint_plot

# save plots with high resolution
ggsave(
  filename = here("figures/plot_kappaEst_ZL2008.jpeg"),
  plot = kappa_plot, width = 6, height = 6
)
ggsave(
  filename = here("figures/plot_pmemEst_ZL2008.jpeg"),
  plot = pMem_plot, width = 6, height = 6
)

ggsave(
  filename = here("figures/plot_jointRes_ZL2008.jpeg"),
  plot = joint_plot, width = 4*2, height = 4
)

## 4.4) Test hypothesis --------------------------------------------------------
# specify hypothesis
hyp_ZL2008 <- c(
  hyp_kappa_1v2 = "kappa1_setsize1 = kappa1_setsize2",
  hyp_kappa_2v3 = "kappa1_setsize2 = kappa1_setsize3",
  hyp_kappa_3v6 = "kappa1_setsize3 = kappa1_setsize6",
  hyp_theta_1v2 = "theta1_setsize1 = theta1_setsize2",
  hyp_theta_2v3 = "theta1_setsize2 = theta1_setsize3",
  hyp_theta_3v6 = "theta1_setsize3 = theta1_setsize6"
)

# test hypothesis
hypothesis(fit_ZL_mixModel,hyp_ZL2008)
