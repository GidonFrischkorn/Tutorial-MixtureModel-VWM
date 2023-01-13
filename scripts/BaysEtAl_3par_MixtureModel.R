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
data_Bays2009 <- read.table(here("data/Bays2009.txt"), header = T) %>% 
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
head(data_Bays2009)

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
                            kappa7 ~ 1,             # uniform  
                            # specify mixing distributions for distinct item categories
                            nlf(theta1 ~ thetat),   # p_mem
                            theta2 ~ 0 + setsize,  # p_intrusion
                            theta3 ~ 0 + setsize,  # p_intrusion
                            theta4 ~ 0 + setsize,  # p_intrusion
                            theta5 ~ 0 + setsize,  # p_intrusion
                            theta6 ~ 0 + setsize,  # p_intrusion
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
                            # thetant ~ 0 + setsize,   # fixed intercept for p_intrusion
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
  # prior(logistic(0, 1), class = b, nlpar = "thetant") +
  prior(constant(-100), class = b, dpar = "theta2", coef = "setsize1") +
  prior(constant(-100), class = b, dpar = "theta3", coef = "setsize1") +
  prior(constant(-100), class = b, dpar = "theta4", coef = "setsize1") +
  prior(constant(-100), class = b, dpar = "theta5", coef = "setsize1") +
  prior(constant(-100), class = b, dpar = "theta6", coef = "setsize1") +
  prior(constant(-100), class = b, dpar = "theta3", coef = "setsize2") +
  prior(constant(-100), class = b, dpar = "theta4", coef = "setsize2") +
  prior(constant(-100), class = b, dpar = "theta5", coef = "setsize2") +
  prior(constant(-100), class = b, dpar = "theta6", coef = "setsize2") +
  prior(constant(-100), class = b, dpar = "theta5", coef = "setsize4") +
  prior(constant(-100), class = b, dpar = "theta6", coef = "setsize4") 



fit_Bays_mixMod <- brm(formula = Bays_mixModel_formula, 
                       data = data_Bays2009, 
                       family = Bays_mixModel, 
                       prior = Bays_mixModel_priors, 
                       iter = 100)


# extract fixed effects
fixef_Bays <- fixef(fit_Bays_mixMod)

kappa_cols <- grepl("kappa",rownames(fixef_Bays))

# extract kappa estimates
kappa_fixedFX <- fixef_Bays[kappa_cols,]

# convert kappa estimates to absolute scale (radians)
kappa_fixedFX <- exp(kappa_fixedFX)

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
  )

# compute pMem
pMem_s1 <- exp(fixef_Bays["thetat_setsize1","Estimate"]) / (
  exp(fixef_Bays["thetat_setsize1","Estimate"]) + exp(0)
)

pMem_s2 <- exp(fixef_Bays["thetat_setsize2","Estimate"]) / (
  exp(fixef_Bays["thetat_setsize2","Estimate"]) +
    1*exp(fixef_Bays["thetant_setsize2","Estimate"]) +
    exp(0) # guessing prob
)

pMem_s4 <- exp(fixef_Bays["thetat_setsize4","Estimate"]) / (
  exp(fixef_Bays["thetat_setsize4","Estimate"]) +
    2*exp(fixef_Bays["thetant_setsize4","Estimate"]) +
    exp(0) # guessing prob
)

pMem_s6 <- exp(fixef_Bays["thetat_setsize6","Estimate"]) / (
  exp(fixef_Bays["thetat_setsize6","Estimate"]) +
    5*exp(fixef_Bays["thetant_setsize6","Estimate"]) +
    exp(0) # guessing prob
)

# compute pSwap
pSwap_s2 <- (1*exp(fixef_Bays["thetant_setsize2","Estimate"])) / (
  exp(fixef_Bays["thetat_setsize2","Estimate"]) +
    1*exp(fixef_Bays["thetant_setsize2","Estimate"]) +
    exp(0) # guessing prob
)

pSwap_s4 <- (2*exp(fixef_Bays["thetant_setsize4","Estimate"])) / (
  exp(fixef_Bays["thetat_setsize4","Estimate"]) +
    2*exp(fixef_Bays["thetant_setsize4","Estimate"]) +
    exp(0) # guessing prob
)

pSwap_s6 <- (5*exp(fixef_Bays["thetant_setsize6","Estimate"])) / (
  exp(fixef_Bays["thetat_setsize6","Estimate"]) +
    5*exp(fixef_Bays["thetant_setsize6","Estimate"]) +
    exp(0) # guessing prob
)

# compute pGuess
pGuess_s1 <- exp(0) / (
  exp(fixef_Bays["thetat_setsize1","Estimate"]) + exp(0)
)

pGuess_s2 <- exp(0) / (
  exp(fixef_Bays["thetat_setsize2","Estimate"]) +
    1*exp(fixef_Bays["thetant_setsize2","Estimate"]) +
    exp(0) # guessing prob
)

pGuess_s4 <- exp(0) / (
  exp(fixef_Bays["thetat_setsize4","Estimate"]) +
    2*exp(fixef_Bays["thetant_setsize4","Estimate"]) +
    exp(0) # guessing prob
)

pGuess_s6 <- exp(0) / (
  exp(fixef_Bays["thetat_setsize6","Estimate"]) +
    5*exp(fixef_Bays["thetant_setsize6","Estimate"]) +
    exp(0) # guessing prob
)

sum(pMem_s1,pGuess_s1)
sum(pMem_s2,pGuess_s2,pSwap_s2)
sum(pMem_s4,pGuess_s4,pSwap_s4)
sum(pMem_s6,pGuess_s6,pSwap_s6)
