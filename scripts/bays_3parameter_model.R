#' This is the tutorial script for setting up the 3-parameter Bays et al (2009) model
#' 
#' 
#' In this script, you will see:
#'  1) how the model is set up using the brms package, 
#'  2) how a simple version of the model is estimates, and 
#'  3) how the model can be evaluated and results extracted and plotted.

rm(list = ls())

library(tidyverse)
library(stats4)
library(circular)
library(lme4)
library(brms)
library(here)

# Set up parallel sampling of mcmc chains
options(mc.cores =  parallel::detectCores())

# specify the number of samples to run for warm up & after warm up
warmup_samples <- 2000
postwarmup_samples <- 2000

# specify the number of chains
nChains <- 4

#' if the number of user defined chains is larger than the number of cores 
#' on the system than estimate as many chains as there are cores on the system
if (nChains >  parallel::detectCores()) {
  nChains <-  parallel::detectCores()
}

mean_se2 <- function (x, mult = 1.96) 
{
  x <- stats::na.omit(x)
  se <- mult * sqrt(stats::var(x)/length(x))
  mean <- mean(x)
  data.frame(y = mean, ymin = mean - se, ymax = mean + se)
}

# function for mixture model likelihood (MLE)
mixtureLL <- function(dat) {
  # transform x from degrees to radians
  rad = dat[c('y','V1','V2','V4','V5')]
  function(p_correct=2, p_other=1, sigma=9) {
    # trasnform laten prob
    p_c = exp(p_correct)/(exp(p_correct)+exp(p_other)+exp(0))
    p_o = exp(p_other)/(exp(p_correct)+exp(p_other)+exp(0))
    p_g = exp(0)/(exp(p_correct)+exp(p_other)+exp(0))
    # transform the normal sd into radians kappa for circular vonmises concentration parameter
    rad_sigma = sigma * pi /180
    kappa = (1/rad_sigma) ** 2
    l_norm <- brms::dvon_mises(rad$y, mu=0, kappa=kappa)
    l_norm1 <- brms::dvon_mises(rad$y, mu=rad$V1, kappa=kappa)
    l_norm2 <- brms::dvon_mises(rad$y, mu=rad$V2, kappa=kappa)
    l_norm4 <- brms::dvon_mises(rad$y, mu=rad$V4, kappa=kappa)
    l_norm5 <- brms::dvon_mises(rad$y, mu=rad$V5, kappa=kappa)
    l_unif <- brms::dvon_mises(rad$y, mu=0, kappa=0)
    likelihood <- p_c*l_norm + p_o/4*(l_norm1+l_norm2+l_norm4+l_norm5) + p_g*l_unif
    -sum(log(likelihood))
  }
}

# fit and return parameter estimates as a mixture of 6 distributions - 1 for correct, 4 for failed bindings, 1 uniform for guessing
fit_mixture3 <- function(dat, init_values=list(p_correct=3, p_other=1, sigma=9)) {
  require(stats4)
  LL_resp <- mixtureLL(dat) 
  # debug(LL_resp)
  fit <- mle(LL_resp, start = init_values)
  coef <- data.frame(t(fit@coef))
  p_correct <- exp(coef$p_correct)/(exp(coef$p_correct)+exp(coef$p_other)+exp(0))
  p_other <- exp(coef$p_other)/(exp(coef$p_correct)+exp(coef$p_other)+exp(0))
  coef$p_correct <- p_correct
  coef$p_other <- p_other
  coef$p_guess <- 1-coef$p_correct-coef$p_other
  coef$negll <- summary(fit)@m2logL
  return(round(coef,3))
}

# load data and transform to radians
dat <- read.csv(here('data/popov_so_reder_exp2and3.csv'))
dat$y <- dat$y * pi /180
dat$V1 <- dat$V1 * pi /180
dat$V2 <- dat$V2 * pi /180
dat$V4 <- dat$V4 * pi /180
dat$V5 <- dat$V5 * pi /180

# fit MLE model
mle_fit <- dat %>% 
  group_by(subject, duration) %>% 
  do({fit_mixture3(.)})

lme4::lmer(p_correct ~ duration + (1|subject), data=mle_fit) %>% summary()

#############################################################################!
# BRMS                                                                   ####
#############################################################################!

# create mixture of von Mises distributions
mix_vonMises1 <- mixture(von_mises(link="identity"),
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
                  # kappa for guessing distribution will be fixed using priors
                  kappa6 ~ 1,         # uniform  
                  # specify mixing distributions for distinct item categories
                  nlf(theta1 ~ thetat),   # p_mem
                  nlf(theta2 ~ thetant),  # p_intrusion
                  nlf(theta3 ~ thetant),  # p_intrusion
                  nlf(theta4 ~ thetant),  # p_intrusion
                  nlf(theta5 ~ thetant),  # p_intrusion
                  # target & guessing distribution will be centered using priors
                  mu1 ~ 1, # fixed intercept constrained using priors
                  mu6 ~ 1, # fixed intercept constrained using priors
                  # center non-target distribution on data specified locations
                  nlf(mu2 ~ V1),           # center non-target
                  nlf(mu3 ~ V2),           # center non-target
                  nlf(mu4 ~ V4),           # center non-target
                  nlf(mu5 ~ V5),           # center non-target
                  # now predict parameters of interest
                  kappa ~ duration + (duration || subject),     # fixed intercept for precision of memory distributions
                  thetat ~ duration + (duration || subject),    # fixed intercept for p_mem
                  thetant ~ duration + (duration || subject),  # fixed intercept for p_intrusion
                  # for brms to process this formula correclty, set non-linear to TRUE
                  nl = TRUE)


# check default priors
get_prior(bf_mixture1, dat, mix_vonMises1)

# constrain priors to identify the model
mix_priors1 <- 
  # first we center the target and guessing distribution to zero
  prior(constant(0), class = Intercept, dpar = "mu1") +
  prior(constant(0), class = Intercept, dpar = "mu6") +
  # next, we set the guessing distribution to be uniform, kappa -> 0
  prior(constant(-100), class = Intercept, dpar = "kappa6")
  # next, we set reasonable priors for the to be estimated distributions
  # prior(normal(5.0, 0.8), class = b, coef = "Intercept", nlpar = "kappa")
  # prior(logistic(0, 1), class = b, coef = "Intercept", nlpar = "theta_t") +
  # prior(logistic(0, 1), class = b, coef = "Intercept", nlpar = "theta_nt") +
  # prior(logistic(0, 1), class = b, coef = "Intercept", nlpar = "theta_g") 


fit1 <- brm(bf_mixture1, dat, mix_vonMises1, mix_priors1, iter = 100, chains=nChains)

exp(fixef(fit1)['kappa_Intercept','Estimate'])

exp(fixef(fit1)['theta1_Intercept','Estimate'])/(exp(fixef(fit1)['theta1_Intercept','Estimate'])+4*exp(fixef(fit1)['thetant_Intercept','Estimate'])+1)
4*exp(fixef(fit1)['thetant_Intercept','Estimate'])/(exp(fixef(fit1)['theta1_Intercept','Estimate'])+4*exp(fixef(fit1)['thetant_Intercept','Estimate'])+1)
1/(exp(fixef(fit1)['theta1_Intercept','Estimate'])+4*exp(fixef(fit1)['thetant_Intercept','Estimate'])+1)