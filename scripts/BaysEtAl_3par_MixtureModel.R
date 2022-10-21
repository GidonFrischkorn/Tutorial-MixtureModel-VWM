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
max_treedepth <- 10

# load data file
data_ZL2008 <- read.table(here("data/Zhang&Luck2008.txt"), header = T) %>% 
  dplyr::mutate(# wrap cases smaller than -pi,
                # or larger than pi around the circle
                RespErr = dplyr::case_when(RespErr < -pi ~ RespErr + 2*pi,
                                           RespErr > pi ~ RespErr - 2*pi,
                                           TRUE ~ RespErr),
                LureIdx1 = case_when(setsize >= 2 ~ 0,
                                     TRUE ~ -100),
                LureIdx2 = case_when(setsize >= 3 ~ 0,
                                     TRUE ~ -100),
                LureIdx3 = case_when(setsize >= 4 ~ 0,
                                     TRUE ~ -100),
                LureIdx4 = case_when(setsize >= 5 ~ 0,
                                     TRUE ~ -100),
                LureIdx5 = case_when(setsize >= 6 ~ 0,
                                     TRUE ~ -100),
                Pos_Lure1 = case_when(is.na(Pos_Lure1) ~ 0,
                                      TRUE ~ Pos_Lure1),
                Pos_Lure2 = case_when(is.na(Pos_Lure2) ~ 0,
                                      TRUE ~ Pos_Lure2),
                Pos_Lure3 = case_when(is.na(Pos_Lure3) ~ 0,
                                      TRUE ~ Pos_Lure3),
                Pos_Lure4 = case_when(is.na(Pos_Lure4) ~ 0,
                                      TRUE ~ Pos_Lure4),
                Pos_Lure5 = case_when(is.na(Pos_Lure5) ~ 0,
                                      TRUE ~ Pos_Lure5),
                setsize = as.factor(setsize),) %>% 
  select(subID,trial,setsize,RespErr,
         Pos_Lure1, Pos_Lure2, Pos_Lure3, Pos_Lure4, Pos_Lure5,
         LureIdx1, LureIdx2, LureIdx3, LureIdx4, LureIdx5)



# have a look at the data and included variables
head(data_ZL2008)

# load custom functions
file.sources = list.files(path = here("functions"), 
                          pattern = "*.R")
sapply(paste(here(),"functions",file.sources,sep = "/"),source,.GlobalEnv)

###############################################################################!
# 1) MLE fit -------------------------------------------------------------------
###############################################################################!

# load data and transform to radians
dat <- read.csv(here('data/popov_so_reder_exp1.csv'))
dat$y  <- dat$y  * pi / 180
dat$V1 <- dat$V1 * pi / 180
dat$V2 <- dat$V2 * pi / 180
dat$V4 <- dat$V4 * pi / 180
dat$V5 <- dat$V5 * pi / 180

# fit MLE model
fit_3par_mixMod(dat)

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
                            kappa7 ~ 1,              # uniform  
                            # specify mixing distributions for distinct item categories
                            nlf(theta1 ~ thetat),   # p_mem
                            nlf(theta2 ~ LureIdx1 + thetant),  # p_intrusion
                            nlf(theta3 ~ LureIdx2 + thetant),  # p_intrusion
                            nlf(theta4 ~ LureIdx3 + thetant),  # p_intrusion
                            nlf(theta5 ~ LureIdx4 + thetant),  # p_intrusion
                            nlf(theta6 ~ LureIdx5 + thetant),  # p_intrusion
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
                            kappa ~ 1,     # fixed intercept for precision of memory distributions
                            thetat ~ 0 + setsize,    # fixed intercept for p_mem
                            thetant ~ 0 + setsize,   # fixed intercept for p_intrusion
                            #thetag ~ 1,    # fixed intercept for p_guess
                            # for brms to process this formula correclty, set non-linear to TRUE
                            nl = TRUE)


# check default priors
get_prior(Bays_mixModel_formula, data_ZL2008, Bays_mixModel)

# constrain priors to identify the model
Bays_mixModel_priors <- 
  # first we center the target and guessing distribution to zero
  prior(constant(0), class = Intercept, dpar = "mu1") +
  prior(constant(0), class = Intercept, dpar = "mu7") +
  # next, we set the guessing distribution to be uniform, kappa -> 0
  prior(constant(-100), class = Intercept, dpar = "kappa7") +
  # next, we set reasonable priors for the to be estimated distributions
  prior(normal(5.0, 0.8), class = b, coef = "Intercept", nlpar = "kappa") +
  prior(logistic(0, 1), class = b, nlpar = "thetat") +
  prior(logistic(0, 1), class = b, nlpar = "thetant") 


fit_Bays_mixMod <- brm(formula = Bays_mixModel_formula, 
                       data = data_ZL2008, 
                       family = Bays_mixModel, 
                       prior = Bays_mixModel_priors, 
                       iter = 2000)

exp(fixef(fit_Bays_mixMod)['kappa_Intercept','Estimate'])

exp(fixef(fit_Bays_mixMod)['thetat_Intercept','Estimate'])/(exp(fixef(fit_Bays_mixMod)['thetat_Intercept','Estimate']) + 4*exp(fixef(fit_Bays_mixMod)['thetant_Intercept','Estimate']) + 1)
4*exp(fixef(fit_Bays_mixMod)['thetant_Intercept','Estimate'])/(exp(fixef(fit_Bays_mixMod)['thetat_Intercept','Estimate']) + 4*exp(fixef(fit_Bays_mixMod)['thetant_Intercept','Estimate'])+1)
1/(exp(fixef(fit_Bays_mixMod)['thetat_Intercept','Estimate'])+ 4 * exp(fixef(fit_Bays_mixMod)['thetant_Intercept','Estimate']) + 1)
