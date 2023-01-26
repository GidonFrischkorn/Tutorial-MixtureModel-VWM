#' This script tests the 3-parameter model of Bays et al (2009) on a synthetic 
#' dataset to see if it can recover the non-target proportion. Setsize varies from
#' 2 to 6 and the probability of non-target responses and guessing increases with
#' setsize
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
library(tidybayes)

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

dat <- data.frame()
pnt <- c(0.2,0.3,0.4)
pmem <- c(0.7, 0.55, 0.4)
ss = c(2,4,6)
for (i in 1:3) {
  DAT <- gen_bays_3p_data(N=5000, pmem=pmem[i], pnt=pnt[i], kappa=10, setsize=ss[i])
  DAT$setsize=ss[i]
  dat <- bind_rows(dat, DAT)
}
dat$setsize <- as.factor(dat$setsize)
dat <- select(dat, setsize,y,nt1_loc, nt2_loc:nt5_loc)

#############################################################################!
# EXPLORE THE DATA                                                       ####
#############################################################################!

head(dat)

# distribution of errors relative to target location (0)
ggplot(dat, aes(y, color=setsize)) +
  geom_density()

# distribution of errors relative to non-target location
nt_errors <- wrap(dat$y-dat[,3:ncol(dat)])
nt_errors$setsize <- dat$setsize
nt_errors <- gather(nt_errors, lureid, y, -setsize)

ggplot(nt_errors, aes(y, fill=setsize)) +
  geom_histogram() +
  facet_wrap(~setsize)

ggplot(nt_errors, aes(y, fill=setsize)) +
  geom_histogram(aes(y=..density..)) +
  facet_wrap(~setsize)

#############################################################################!
# ADD indices necessary for brms                                         ####
#############################################################################!

dat<- mutate(dat,
  setsize  = as.numeric(as.character(setsize)),
  LureIdx1 = case_when(setsize >= 2 ~ 1,
                       TRUE ~ 0),
  LureIdx2 = case_when(setsize >= 3 ~ 1,
                       TRUE ~ 0),
  LureIdx3 = case_when(setsize >= 4 ~ 1,
                       TRUE ~ 0),
  LureIdx4 = case_when(setsize >= 5 ~ 1,
                       TRUE ~ 0),
  LureIdx5 = case_when(setsize >= 6 ~ 1,
                       TRUE ~ 0),
  inv_ss = 1/(setsize-1),
  setsize = as.factor(setsize))

dat[is.na(dat)] <- 0

#############################################################################!
# BRMS                                                                   ####
#############################################################################!


if (!file.exists(here('output/fit_3pmixt_recovery_var_ss_logconst_larger_dataset.RData'))) {
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
                    nlf(theta2 ~ LureIdx1*(thetant + log(inv_ss)) + (1-LureIdx1)*(-100)),   # p_intrusion
                    nlf(theta3 ~ LureIdx2*(thetant + log(inv_ss)) + (1-LureIdx2)*(-100)),   # p_intrusion
                    nlf(theta4 ~ LureIdx3*(thetant + log(inv_ss)) + (1-LureIdx3)*(-100)),   # p_intrusion
                    nlf(theta5 ~ LureIdx4*(thetant + log(inv_ss)) + (1-LureIdx4)*(-100)),   # p_intrusion
                    nlf(theta6 ~ LureIdx5*(thetant + log(inv_ss)) + (1-LureIdx5)*(-100)),   # p_intrusion
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
                    thetat ~ 0+setsize,              # fixed intercept for p_mem
                    thetant ~ 0+setsize,             # fixed intercept for p_intrusion
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
    prior(logistic(0, 1), class = b, nlpar = "thetat") +
    prior(logistic(0, 1), class = b, nlpar = "thetant")
  
  
  fit <- brm(bf_mixture1, dat, mix_vonMises1, mix_priors1, iter = warmup_samples+postwarmup_samples, chains=nChains)
  
  save(list=ls(), file=here('output/fit_3pmixt_recovery_var_ss_logconst_larger_dataset.RData'))
} else {
  load(here('output/fit_3pmixt_recovery_var_ss_logconst_larger_dataset.RData'))
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
theta <- theta[,1]

# the model recovers both kappa and the proportions well
exp(kappa)

pmem_est <- c(0,0,0)
pnt_est <- c(0,0,0)

pmem_est[1] <- exp(theta[1])/sum(exp(theta[1])+exp(theta[4])+exp(0))
pnt_est[1] <- exp(theta[4])/sum(exp(theta[1])+exp(theta[4])+exp(0))
exp(0)/sum(exp(theta[1])+exp(theta[4])+exp(0))


pmem_est[2] <- exp(theta[2])/c(exp(theta[2])+exp(theta[5])+exp(0))
pnt_est[2] <- exp(theta[5])/c(exp(theta[2])+exp(theta[5])+exp(0))
exp(0)/c(exp(theta[2])+exp(theta[5])+exp(0))

pmem_est[3] <- exp(theta[3])/c(exp(theta[3])+exp(theta[6])+exp(0))
pnt_est[3] <- exp(theta[6])/c(exp(theta[3])+exp(theta[6])+exp(0))
exp(0)/c(exp(theta[3])+exp(theta[6])+exp(0))

pmem_est
pnt_est
# true
#pnt <- c(0.2,0.3,0.4)
#pmem <- c(0.7, 0.55, 0.4)



nt_errors <- wrap(dat$y-dat[,3:ncol(dat)])
nt_errors$setsize <- dat$setsize
nt_errors <- gather(nt_errors, lureid, y, -setsize)

ggplot(nt_errors, aes(y, fill=setsize)) +
  geom_histogram() +
  facet_wrap(~setsize)

ggplot(nt_errors, aes(y, fill=setsize)) +
  geom_histogram(aes(y=..density..)) +
  facet_wrap(~setsize)

pp <- pp_check(fit)
pp <- pp$data
pp <- arrange(pp, rep_id, y_id)

