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
pacman::p_load(here, brms, tidyverse, tidybayes, patchwork, gghalves, bmm)

# load function to clean up plots
source(here("functions","clean_plot.R"))

# load missing output files
source(here("scripts","LoadResultsFiles.R"))

# Set up parallel sampling of mcmc chains
options(mc.cores =  parallel::detectCores())

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

# generate random imm parameters to simulate data from
# simulation settings
nsubj = 30
ntrials = 200
setsize = 5

# random parameters for each subject
parms <- data.frame(
  c = rnorm(n = nsubj, mean = 2, sd = 0.5),
  a = rnorm(n = nsubj, mean = 0.5, sd = 0.15),
  n = 0,
  s = pmax(0,rnorm(n = nsubj, mean = 2, sd = 0.5)),
  kappa = pmax(0,rnorm(nsubj, mean = 5, sd = 1.5))
)

# simulate data from the IMM parameters
simData <- gen_imm_data(parms = parms, ntrial = ntrials, setsize = setsize)

# add indices for setsize
simData <- simData %>% 
  dplyr::mutate(
    LureIdx1 = case_when(setsize >= 2 ~ 1,
                         TRUE ~ 0),
    LureIdx2 = case_when(setsize >= 3 ~ 1,
                         TRUE ~ 0),
    LureIdx3 = case_when(setsize >= 4 ~ 1,
                         TRUE ~ 0),
    LureIdx4 = case_when(setsize >= 5 ~ 1,
                         TRUE ~ 0),
    inv_SS = 1/(setsize - 1),
    setsize = as.factor(setsize)) %>% 
  select(ID,trial,setsize,respErr,inv_SS,
         Item1_rel, Item2_rel, Item3_rel, Item4_rel, Item5_rel,
         spaD1, spaD2, spaD3, spaD4, spaD5,
         LureIdx1,LureIdx2,LureIdx3,LureIdx4)

###############################################################################!
# 2) BRMS fit ------------------------------------------------------------------
###############################################################################!

# create mixture of von Mises distributions
IMM_mixModel <- mixture(von_mises(link = "identity"),
                        von_mises(link = "identity"),
                        von_mises(link = "identity"),
                        von_mises(link = "identity"),
                        von_mises(link = "identity"),
                        von_mises(link = "identity"),
                        order = "none")

# set up mixture model
IMM_mixModel_formula <- bf(respErr ~ 1,
                           # fix kappa over memory distributions
                           nlf(kappa1 ~ kappa),     # target distribution
                           nlf(kappa2 ~ kappa),     # non-target
                           nlf(kappa3 ~ kappa),     # non-target
                           nlf(kappa4 ~ kappa),     # non-target
                           nlf(kappa5 ~ kappa),     # non-target
                           # kappa for guessing distribution will be fixed using priors
                           kappa6 ~ 1,             # uniform  
                           # specify mixing distributions for distinct item categories
                           nlf(theta1 ~ exp(-s*spaD1)*c + a),   # p_mem
                           nlf(theta2 ~ LureIdx1*(exp(-s*spaD2)*c + a) + (1-LureIdx1)*(-100)),  # p_intrusion
                           nlf(theta3 ~ LureIdx2*(exp(-s*spaD3)*c + a) + (1-LureIdx2)*(-100)),  # p_intrusion
                           nlf(theta4 ~ LureIdx3*(exp(-s*spaD4)*c + a) + (1-LureIdx3)*(-100)),  # p_intrusion
                           nlf(theta5 ~ LureIdx4*(exp(-s*spaD5)*c + a) + (1-LureIdx4)*(-100)),  # p_intrusion
                           # target & guessing distribution will be centered using priors
                           mu1 ~ 1, # fixed intercept constrained using priors
                           mu6 ~ 1, # fixed intercept constrained using priors
                           # center non-target distribution on data specified locations
                           nlf(mu2 ~ Item2_rel),           # center non-target
                           nlf(mu3 ~ Item3_rel),           # center non-target
                           nlf(mu4 ~ Item4_rel),           # center non-target
                           nlf(mu5 ~ Item5_rel),           # center non-target
                           nlf(s ~ exp(logS)),      
                           # now predict parameters of interest
                           kappa ~ 1 + (1 || ID),  # fixed intercept & random slope: precision of memory distributions
                           logS ~ 1 + (1 || ID),   # fixed intercept & random slope: spatial gradient (on logarithmic scale)
                           c ~ 1 + (1 || ID),      # fixed intercept & random slope: context activation
                           a ~ 1 + (1 || ID),      # fixed intercept & random slope: general activation
                           # for brms to process this formula correctly, set non-linear to TRUE
                           nl = TRUE)

# check default priors
get_prior(IMM_mixModel_formula, simData, IMM_mixModel)
 
# constrain priors to identify the model
IMM_priors <- 
  # first we center the target and guessing distribution to zero
  prior(constant(0), class = Intercept, dpar = "mu1") + 
  prior(constant(0), class = Intercept, dpar = "mu6") +
  # next, we set the guessing distribution to be uniform, kappa -> 0
  prior(constant(-100), class = Intercept, dpar = "kappa6") +
  # next, we set reasonable priors for the to be estimated distributions
  prior(normal(5.0, 0.8), class = b, nlpar = "kappa") +
  prior(normal(0, 1), class = b, nlpar = "logS", ) +
  prior(normal(2, 1), class = b, nlpar = "c") +
  prior(normal(0.5, 0.2), class = b, nlpar = "a") 



fit_IMM_mixMod <- brm(formula = IMM_mixModel_formula, 
                      data = simData,
                      family = IMM_mixModel, 
                      prior = IMM_priors, 
                      iter = 1000 + 1000, 
                      chains = nChains,
                      save_pars = save_pars(all = T))

summary(fit_IMM_mixMod)
apply(parms,2,mean)

randFX <- ranef(fit_IMM_mixMod)

cor(randFX$ID[,"Estimate","c_Intercept"],parms$c)
cor(randFX$ID[,"Estimate","a_Intercept"],parms$a)
cor(randFX$ID[,"Estimate","logS_Intercept"],log(parms$s))
cor(randFX$ID[,"Estimate","kappa_Intercept"],parms$kappa)
