#' This is the first tutorial script for setting up the Zhang & Luck (2008) mixture model
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
library(ggplot2)
library(magrittr)
library(brmsmargins)

# Set up parallel sampling of mcmc chains
options(mc.cores =  parallel::detectCores())

# specify the number of samples to run for warm up & after warm up
warmup_samples <- 1000
postwarmup_samples <- 1000

# specify the number of chains
nChains <- 10

#' if the number of user defined chains is larger than the number of cores 
#' on the system than estimate as many chains as there are cores on the system
if (nChains >  parallel::detectCores()) {
  nChains <-  parallel::detectCores()
}

# set brms controls to prevent divergent transitions & improve convergence
adapt_delta <- .95
max_treedepth <- 12

# load data file
data_ZL2008 <- read.table(paste(getwd(),"data","Zhang&Luck2008.txt", sep = "/"), header = T) %>% 
  dplyr::mutate(setsize = as.factor(setsize),
                # wrap cases smaller than -pi,
                # or larger than pi aournd the circle
                RespErr = dplyr::case_when(RespErr < -pi ~ RespErr + 2*pi,
                                           RespErr > pi ~ RespErr - 2*pi,
                                           TRUE ~ RespErr))

# have a look at the data and included variables
head(data_ZL2008)

# 1) Model Setup ---------------------------------------------------------------

#' First we specify a mixture of von Mises distributions.
#' The first distribution is the memory distribution and
#' the second distribution is the guessing
ZL_mixModel <- mixture(von_mises,von_mises, order = "none")

#' Then, we set up the formula for the mixture model.
#' Although, we do not want to estimate the mean of the two von Mises distributions,
#' we have to initialize the formula to specify the dependent variable.
#' Using priors (see step 3), we will constrain the means of both von Mises distributions
#' to zero. Additionally, we will use priors to fix the precision (kappa) of the
#' second von Mises to be zero (at least practically zero). 
ZL_mixModel_formula <- bf(RespErr ~ 1,    # Initializing the dependent variable
                          # estimating fixed intercept & random intercept for kappa of the first von Mises
                          kappa1 ~ 1 + (1 || subID), 
                          kappa2 ~ 1,
                          # estimating fixed intercept & random intercept for the mixing proportion 
                          # of the first vonMises (i.e., p_mem)
                          theta1 ~ 1 + (1 || subID)) 

# constrain parameters using priors
ZL_mixModel_priors <- 
  # fix mean of the first von Mises to zero
  prior(constant(0), class = Intercept, dpar = "mu1") +
  # fix mean of the second von Mises to zero
  prior(constant(0), class = Intercept, dpar = "mu2") +
  # fix kappa of the second von Mises to (alomst) zero
  prior(constant(-100), class = Intercept, dpar = "kappa2")

# 2) Model estimation ----------------------------------------------------------

# fit mixture model if there is not already a results file stored
if (!file.exists(paste(getwd(),"output","fit_mixtureMod_ZL2008.RData", sep = "/"))) {
  #' To reproduce the results from Zhang & Luck (2008), we include setsize
  #' as a within subject predictor.
  #' Additionally we specify the formula so that we directly get the estimates
  #' for each level of setsize: 0 + setsize
  #' Finally, we implement the hierarchical structure by allowing that the setsize
  #' effect varies over each subject: (0 + setsize || subID)
  #' This is done for both kappa1, the precision of the memory distribution,
  #' and theta1, the mixing distribution, essentially estimating pMem.
  ZL_mixModel_formula <- bf(RespErr ~ 1,    # Initializing the dependent variable
                            # estimating fixed intercept & random intercept for kappa of the first von Mises
                            kappa1 ~ 0 + setsize + (0 + setsize || subID), 
                            kappa2 ~ 1,
                            # estimating fixed intercept & random intercept for the mixing proportion 
                            # of the first vonMises (i.e., p_mem)
                            theta1 ~ 0 + setsize + (0 + setsize || subID))
  
  
  # using this adapted model formula we can now fit the mixture model using brms
  fit_mixtureMod <- brm(
    # include model information
    formula = ZL_mixModel_formula, # specify formula for mixture model
    data    = data_ZL2008, # specify data used to estimate the mixture model
    family  = ZL_mixModel, # call the defined mixture family
    prior   = ZL_mixModel_priors, # use the used defined priors,
    
    # add brms settings
    warmup = warmup_samples,
    iter = warmup_samples + postwarmup_samples, 
    chains = nChains,
    control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
  )
  
  # save results into file
  save(fit_mixtureMod, 
       file = paste(getwd(),"output","fit_mixtureMod_ZL2008.RData", sep = "/"))
  
} else {
  # load results file
  load(file = paste(getwd(),"output","fit_mixtureMod_ZL2008.RData", sep = "/"))
}

# 3) Model evaluation ----------------------------------------------------------

# plot the posterior predicitive check to evaluate overall model fit
pp_check(fit_mixtureMod)

# print results summary
summary(fit_mixtureMod)



# extract the fixed effects from the model
fixedEff <- fixef(fit_mixtureMod)

# determine the rows that contain the relevant parameter estimates
theta_cols <- grepl("theta",rownames(fixedEff))
kappa_cols <- grepl("kappa1",rownames(fixedEff))

# extract kappa estimates
kappa_fixedFX <- fixedEff[kappa_cols,]

# prepare kappa estimates for plotting
df_kappa_plot <- as.data.frame(exp(kappa_fixedFX)) %>% 
  dplyr::mutate(sd_rad = sqrt(1/Estimate),
          sd_rad_UL = sqrt(1/Q2.5),  # lower precision is higher s.d. 
          sd_rad_LL = sqrt(1/Q97.5), # higher precision is lower s.d.
          setsize = levels(data$setsize))

# plot s.d. estimates over setsizes
ggplot(data = df_kappa_plot,
       aes(x = setsize, y = sd_rad, ymin = sd_rad_LL, ymax = sd_rad_UL)) +
  geom_pointrange() +
  labs(x = "Memory Set Size", y = "s.d. (in Radians)") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank())

# extract theta estimates
theta_fixedFX <- fixedEff[theta_cols,]

# convert theta estimates into pMem estimates
p_Mem_fixedFX <- gtools::inv.logit(theta_fixedFX)

# prepare pMem estimates for ploting
df_pMem_plot <- as.data.frame(p_Mem_fixedFX) %>% 
  dplyr::mutate(setsize = levels(data$setsize))

# plot pMem estimates across setsize
ggplot(data = df_pMem_plot,
       aes(x = setsize, y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
  geom_pointrange() + 
  labs(x = "Memory Set Size", y = "pMem") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank())
