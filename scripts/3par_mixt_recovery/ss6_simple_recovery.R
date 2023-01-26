#' This script tests the 3-parameter model of Bays et al (2009) on a synthetic 
#' dataset to see if it can recover the non-target proportion. Setsize is fixed
#' to 6
#' 
#' @author: Ven Popov
#' 

#############################################################################!
# Preliminaries                                                          ####
#############################################################################!

rm(list = ls())

library(tidyverse)
library(stats4)
library(circular)
library(lme4)
library(brms)
library(here)

source(here('functions/wrap.R'))
source(here('functions/gen_bays_3p_data.R'))

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

#############################################################################!
# Generate synthetic data                                                ####
#############################################################################!

dat <- gen_bays_3p_data(N=2000, pmem=0.6, pnt=0.3, kappa=10, setsize=6)

head(dat)

# distribution of errors relative to target location (0)
hist(wrap(dat$y), breaks=60)

# distribution of errors relative to non-target location
nt_errors <- as.matrix(wrap(dat$y-dat[,2:ncol(dat)]))
hist(nt_errors, breaks=60)

#############################################################################!
# BRMS                                                                   ####
#############################################################################!


if (!file.exists(here('output/fit_3pmixt_recovery_ss6.RData'))) {
  # create mixture of von Mises distributions
  mix_vonMises1 <- mixture(von_mises(link="identity"),
                           von_mises(link="identity"),
                           von_mises(link="identity"),
                           von_mises(link="identity"),
                           von_mises(link="identity"),
                           von_mises(link="identity"),
                           von_mises(link="identity"),
                           order = "none")
  
  # set up mixture model
  bf_mixture1 <- bf(y ~ 1,
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
                    nlf(theta1 ~ thetat),    # p_mem
                    nlf(theta2 ~ thetant),   # p_intrusion
                    nlf(theta3 ~ thetant),   # p_intrusion
                    nlf(theta4 ~ thetant),   # p_intrusion
                    nlf(theta5 ~ thetant),   # p_intrusion
                    nlf(theta6 ~ thetant),   # p_intrusion
                    # target & guessing distribution will be centered using priors
                    mu1 ~ 1, # fixed intercept constrained using priors
                    mu7 ~ 1, # fixed intercept constrained using priors
                    # center non-target distribution on data specified locations
                    nlf(mu2 ~ nt1_loc),      # center non-target
                    nlf(mu3 ~ nt2_loc),      # center non-target
                    nlf(mu4 ~ nt3_loc),      # center non-target
                    nlf(mu5 ~ nt4_loc),      # center non-target
                    nlf(mu6 ~ nt5_loc),      # center non-target
                    # now predict parameters of interest
                    kappa ~ 1,               # fixed intercept for precision of memory distributions
                    thetat ~ 1,              # fixed intercept for p_mem
                    thetant ~ 1,             # fixed intercept for p_intrusion
                    # for brms to process this formula correclty, set non-linear to TRUE
                    nl = TRUE)
  
  
  # check default priors
  get_prior(bf_mixture1, dat, mix_vonMises1)
  
  # constrain priors to identify the model
  mix_priors1 <- 
    # first we center the target and guessing distribution to zero
    prior(constant(0), class = Intercept, dpar = "mu1") +
    prior(constant(0), class = Intercept, dpar = "mu7") +
    # next, we set the guessing distribution to be uniform, kappa -> 0
    prior(constant(-100), class = Intercept, dpar = "kappa7") +
    # next, we set reasonable priors for the to be estimated distributions
    prior(normal(5.0, 0.8), class = b, coef = "Intercept", nlpar = "kappa") +
    prior(logistic(0, 1), class = b, coef = "Intercept", nlpar = "thetat") +
    prior(logistic(0, 1), class = b, coef = "Intercept", nlpar = "thetant")
  
  
  fit <- brm(bf_mixture1, dat, mix_vonMises1, mix_priors1, iter = warmup_samples+postwarmup_samples, chains=nChains)
  
  save(list=ls(), file=here('output/fit_3pmixt_recovery_ss6.RData'))
} else {
  load(here('output/fit_3pmixt_recovery_ss6.RData'))
}

#############################################################################!
# RESULTS                                                                ####
#############################################################################!

# view estimates
summary(fit)

# extract parameters
fixedF <- fixef(fit)
kappa <- fixedF['kappa_Intercept',]
theta <- fixedF[grepl('theta',row.names(fixedF)),]
theta <- c(theta[,1],0)

# the model recovers both kappa and the proportions well
exp(kappa)
c(1,5,1)*exp(theta)/sum(exp(theta)*c(1,5,1))
