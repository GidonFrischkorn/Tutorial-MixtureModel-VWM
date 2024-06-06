 #' This is the tutorial script for setting up the Zhang & Luck (2008) mixture model
#' for visual working memory tasks that use continuous report recall procedures.
#' 
#' In this script, you will see:
#'  1) how the model is set up using the bmm package, 
#'  2) how a simple version of the model is estimates, and 
#'  3) how the model can be evaluated and results extracted and plotted.

###############################################################################!
# 0) R Setup: Packages & Data --------------------------------------------------
###############################################################################!

# start fresh
rm(list = ls()) # clean up work space
graphics.off()  # switch off graphics device

# load required packages
pacman::p_load(here, brms, tidyverse, tidybayes, patchwork, gghalves)
pacman::p_load_gh("venpopov/bmm")

# load missing output files
source(here("scripts","LoadResultsFiles.R"))

# load function to clean up plots
source(here("functions","clean_plot.R"))

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
adapt_delta <- .99
max_treedepth <- 12

###############################################################################!
# 1) Read: Raw Data ------------------------------------------------------------
###############################################################################!

# load data file
data_ZL2008 <- zhang_luck_2008

# have a look at the data and included variables
head(data_ZL2008)

###############################################################################!
# 2) Model Setup ---------------------------------------------------------------
###############################################################################!

#' First, we set up the bmmodel object to specify that we want to fit a two-parameter
#' mixture model and connect the relevant variables from our data to the model
ZL_model <- mixture2p(resp_error = "response_error")

#' To reproduce the results from Zhang & Luck (2008), we include setsize
#' as a within subject predictor.
#' Additionally we specify the formula so that we directly get the estimates
#' for each level of setsize: 0 + setsize
#' Finally, we implement the hierarchical structure by allowing that the setsize
#' effect varies over each subject: (0 + setsize || subID)
#' This is done for both kappa, the precision of the memory distribution,
#' and thetat, the mixing distribution for target responses, essentially estimating pMem.
ZL_formula <- bmf(
  # estimating fixed intercept & random intercept for kappa of the first von Mises
  kappa ~ 0 + setsize + (0 + setsize || subID), 
  # estimating fixed intercept & random intercept for the mixing proportion 
  # for the memory/target distribution (i.e., p_mem)
  thetat ~ 0 + setsize + (0 + setsize || subID)
)

# We can access the default priors generated for this model using the default_prior function
default_prior(ZL_formula, data = data_ZL2008, model = ZL_model)

###############################################################################!
# 3) Model estimation ----------------------------------------------------------
###############################################################################!

# using the bmm function, we pass the bmm formula for the mixture model and
# the specified bmmodel. All other constraints are taken care of within this function
# additionally any further arguments for brms can be passed as well
ZL_fit <- bmm(
  formula = ZL_formula, # specify formula for mixture model
  data    = data_ZL2008,   # specify data used to estimate the mixture model
  model = ZL_model, # select the two-parameter model for fitting
  
  # save all potentially relevant information
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
  
  backend = "cmdstanr",
  refresh = 100,
  file = here("output/fit_E1_ZL2008_bmm")
)


###############################################################################!
# 4) Model evaluation ----------------------------------------------------------
###############################################################################!

## 4.1) fit & summary ----------------------------------------------------------

round(max(rhat(ZL_fit), na.rm = T),3)
hist(neff_ratio(ZL_fit))
nuts_params(ZL_fit)
hist(log_posterior(ZL_fit)$Value)

# plot the posterior predictive check to evaluate overall model fit
brms::pp_check(ZL_fit, group = "setsize", type = "dens_overlay_grouped")

# quick plot of the conditional effects on the two parameters
conditional_effects(ZL_fit, effects = "setsize", dpar = "kappa1")
conditional_effects(ZL_fit, effects = "setsize", nlpar = "kappa")

conditional_effects(ZL_fit, effects = "setsize", nlpar = "thetat")
conditional_effects(ZL_fit, effects = "setsize", dpar = "theta1")

# print results summary
summary(ZL_fit)

## 4.2) extract parameter estimates --------------------------------------------

# extract the fixed effects from the model
fixedEff <- fixef(ZL_fit)

# determine the rows that contain the relevant parameter estimates
theta_cols <- startsWith(rownames(fixedEff),"thetat_")
kappa_cols <- startsWith(rownames(fixedEff),"kappa_")

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
fixedFX_draws <- ZL_fit %>% 
  tidy_draws() %>%
  select(starts_with("b_"),.chain,.iteration,.draw) %>% 
  pivot_longer(cols = starts_with("b_"),
               names_to = "modelPar",
               values_to = "postSample") %>% 
  mutate(par = str_split_i(modelPar,"_",2),
         cond = str_split_i(modelPar,"_",3)) %>% 
  select(-modelPar) %>% 
  filter(par == "kappa" | par == "thetat") %>% 
  mutate(postSample_abs = case_when(par == "kappa" ~ (sqrt(1/exp(postSample))/pi) * 180,
                                    par == "thetat" ~ inv_logit_scaled(postSample)),
         cond = str_remove_all(cond,"setsize"))

results_ZL2008 <- data.frame(
  cond = as.character(c(1,2,3,6)),
  pMem = c(.99,.95,.83,.38),
  sdMen = c(13.9,19.4,21.9,22.3)
)

# plot kappa results
kappa_plot <- ggplot(data = fixedFX_draws %>% filter(par == "kappa"),
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
pMem_plot <- ggplot(data = fixedFX_draws %>% filter(par == "thetat"),
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
  hyp_kappa_1v2 = "kappa_setsize1 = kappa_setsize2",
  hyp_kappa_2v3 = "kappa_setsize2 = kappa_setsize3",
  hyp_kappa_3v6 = "kappa_setsize3 = kappa_setsize6",
  hyp_theta_1v2 = "thetat_setsize1 = thetat_setsize2",
  hyp_theta_2v3 = "thetat_setsize2 = thetat_setsize3",
  hyp_theta_3v6 = "thetat_setsize3 = thetat_setsize6"
)

# test hypothesis
hypothesis(ZL_fit,hyp_ZL2008)
plot(hypothesis(ZL_fit,hyp_ZL2008))

# 5) Extract brms info from bmm object -----------------
brmsformula <- ZL_fit$formula
brmsfamily <- ZL_fit$family
brmsdata <- ZL_fit$data

brmsformula
brmsfamily
head(brmsdata)
