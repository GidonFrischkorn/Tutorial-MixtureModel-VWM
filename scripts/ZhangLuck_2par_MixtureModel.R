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
library(brms)       # for estimating the mixture model
library(tidyverse)  # for data wrangling & ploting
library(here)       # use of relative file paths

# Set up parallel sampling of mcmc chains
options(mc.cores =  parallel::detectCores())

# specify the number of samples to run for warm up & after warm up
warmup_samples <- 2000
postwarmup_samples <- 2000

# specify the number of chains
nChains <- 5

#' if the number of user defined chains is larger than the number of cores 
#' on the system than estimate as many chains as there are cores on the system
if (nChains >  parallel::detectCores()) {
  nChains <-  parallel::detectCores()
}

# set brms controls to prevent divergent transitions & improve convergence
adapt_delta <- .95
max_treedepth <- 10

# load data file
data_ZL2008 <- read.table(here("data/Zhang&Luck2008.txt"), header = T) %>% 
  dplyr::mutate(setsize = as.factor(setsize),
                # wrap cases smaller than -pi,
                # or larger than pi around the circle
                RespErr = dplyr::case_when(RespErr < -pi ~ RespErr + 2*pi,
                                           RespErr > pi ~ RespErr - 2*pi,
                                           TRUE ~ RespErr))

# have a look at the data and included variables
head(data_ZL2008)

###############################################################################!
# 1) Model Setup ---------------------------------------------------------------
###############################################################################!

#' First we specify a mixture of von Mises distributions.
#' The first distribution is the memory distribution and
#' the second distribution is the guessing
ZL_mixFamily <- mixture(von_mises,von_mises, order =  "none")

#' Then, we set up the formula for the mixture model.
#' Although, we do not want to estimate the mean of the two von Mises distributions,
#' we have to initialize the formula to specify the dependent variable.
#' Using priors (see step 3), we will constrain the means of both von Mises distributions
#' to zero. Additionally, we will use priors to fix the precision (kappa) of the
#' second von Mises to be zero (at least practically zero). 
ZL_mixFormula <- bf(RespErr ~ 1,    # Initializing the dependent variable
                    mu2 ~ 1,
                    # estimating fixed intercept & random intercept for kappa of the first von Mises
                    kappa1 ~ 1,
                    kappa2 ~ 1,
                    # estimating fixed intercept & random intercept for the mixing proportion 
                    # of the first vonMises (i.e., p_mem)
                    theta1 ~ 1) 

get_prior(formula = ZL_mixFormula,
          family = ZL_mixFamily,
          data = data_ZL2008)

# constrain parameters using priors
ZL_mixPriors <- 
  # fix mean of the first von Mises to zero
  prior(constant(0), class = Intercept, dpar = "mu1") +
  # fix mean of the second von Mises to zero
  prior(constant(0), class = Intercept, dpar = "mu2") +
  # fix kappa of the second von Mises to (alomst) zero
  prior(constant(-100), class = Intercept, dpar = "kappa2")

###############################################################################!
# 2) Model estimation ----------------------------------------------------------
###############################################################################!

# fit mixture model if there is not already a results file stored
if (!file.exists(here("output/fit_ZL2008_mixModel.RData"))) {
  #' To reproduce the results from Zhang & Luck (2008), we include setsize
  #' as a within subject predictor.
  #' Additionally we specify the formula so that we directly get the estimates
  #' for each level of setsize: 0 + setsize
  #' Finally, we implement the hierarchical structure by allowing that the setsize
  #' effect varies over each subject: (0 + setsize || subID)
  #' This is done for both kappa1, the precision of the memory distribution,
  #' and theta1, the mixing distribution, essentially estimating pMem.
  ZL_mixFormula <- bf(RespErr ~ 1,    # Initializing the dependent variable
                      mu2 ~ 1,
                      # estimating fixed intercept & random intercept for kappa of the first von Mises
                      kappa1 ~ 0 + setsize + (0 + setsize || subID), 
                      kappa2 ~ 1,
                      # estimating fixed intercept & random intercept for the mixing proportion 
                      # of the first vonMises (i.e., p_mem)
                      theta1 ~ 0 + setsize + (0 + setsize || subID))
  
  
  # using this adapted model formula we can now fit the mixture model using brms
  fit_ZL_mixModel <- brm(
    # include model information
    formula = ZL_mixFormula, # specify formula for mixture model
    data    = data_ZL2008, # specify data used to estimate the mixture model
    family  = ZL_mixFamily, # call the defined mixture family
    prior   = ZL_mixPriors, # use the used defined priors,
    
    # add brms settings
    warmup = warmup_samples,
    iter = warmup_samples + postwarmup_samples, 
    chains = nChains,
    control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
  )
  
  # save results into file
  save(fit_ZL_mixModel, 
       file = here("output/fit_ZL2008_mixModel.RData"))
  
} else {
  # load results file
  load(file = here("output/fit_ZL2008_mixModel.RData"))
}

###############################################################################!
# 3) Model evaluation ----------------------------------------------------------
###############################################################################!

## 3.1) fit & summary ----------------------------------------------------------
# plot the posterior predictive check to evaluate overall model fit
pp_check(fit_ZL_mixModel)

# print results summary
summary(fit_ZL_mixModel)

## 3.2) extract parameter estimates --------------------------------------------

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

## 3.3) plot parameter estimates -----------------------------------------------

# define defaults for clean ggplots
clean_plot <- theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    axis.line.x = element_line(color = 'black'),
                    axis.line.y = element_line(color = 'black'),
                    legend.key = element_rect(fill = 'white'),
                    text = element_text(size = 15),
                    line = element_line(size = 1),
                    axis.ticks = element_line(size = 1))

# prepare kappa estimates for plotting
df_kappa_plot <- as.data.frame(kappa_fixedFX) %>% 
  dplyr::mutate(
    # convert kappa to the standard deviation in radians
    sd_rad = sqrt(1/Estimate),
    sd_rad_UL = sqrt(1/Q2.5),  # lower precision is higher s.d. 
    sd_rad_LL = sqrt(1/Q97.5), # higher precision is lower s.d.
    # convert standard deviation in radians to degrees
    sd_deg = sd_rad / pi * 180,
    sd_deg_UL = sd_rad_UL / pi * 180,
    sd_deg_LL = sd_rad_LL / pi * 180,
    # add set size variable
    setsize = levels(data_ZL2008$setsize),
    ZL_results = c(13.9,19.4,21.9,22.3)
  )

# plot s.d. estimates over set sizes
kappa_plot <- ggplot(data = df_kappa_plot,
                     aes(x = setsize, y = sd_deg, ymin = sd_deg_LL, ymax = sd_deg_UL)) +
  geom_pointrange() +
  geom_point(aes(y = ZL_results), shape = 23, color = "darkred", fill = "darkred", 
             position = position_nudge(x = 0.05)) +
  labs(x = "Set size", y = "s.d.") +
  coord_cartesian(ylim = c(0,40)) +
  clean_plot

# show plot for precision results
kappa_plot

# prepare pMem estimates for ploting
df_pMem_plot <- as.data.frame(p_Mem_fixedFX) %>% 
  dplyr::mutate(setsize = levels(data_ZL2008$setsize),
                ZL_results = c(.99, .95, .83,.38))

# plot pMem estimates across setsize
pMem_plot <- ggplot(data = df_pMem_plot,
                    aes(x = setsize, y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
  geom_pointrange() + 
  geom_point(aes(y = ZL_results), shape = 23, color = "darkred", fill = "darkred", 
             position = position_nudge(x = 0.05)) +
  labs(x = "Set Size", y = expression(P[mem])) +
  coord_cartesian(ylim = c(0,1)) +
  clean_plot

# show plot for pMem results
pMem_plot

# export plots

# the patchwork pacakge allows to easily joint plots
library(patchwork)
joint_plot <- (pMem_plot | kappa_plot)

# show joint plot
joint_plot

# save plots with high resolution
ggsave(
  filename = here("figures/plot_kappaEst_ZL2008.jpeg"),
  plot = kappa_plot, width = 4, height = 4
)
ggsave(
  filename = here("figures/plot_pmemEst_ZL2008.jpeg"),
  plot = pMem_plot, width = 4, height = 4
)

ggsave(
  filename = here("figures/plot_jointRes_ZL2008.jpeg"),
  plot = joint_plot, width = 4*2, height = 4
)
