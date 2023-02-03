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

# read in data for Experiment 2 from Oberauer & Lin (2017)
df_OberauerLin2017_E1 <- read.table(here("data","OberauerLin2017_IM","colorwheel9.dat"))
colnames(df_OberauerLin2017_E1) <- c(
  "ID","Session","Trial","TrialAlt","SetSize",
  "Item1_Col","Item1_Pos","Item2_Col","Item2_Pos","Item3_Col","Item3_Pos","Item4_Col","Item4_Pos",
  "Item5_Col","Item5_Pos","Item6_Col","Item6_Pos","Item7_Col","Item7_Pos","Item8_Col","Item8_Pos",
  "Response"
)

df_OberauerLin2017_E1 <- df_OberauerLin2017_E1 %>% 
  mutate(deviation = (Response - Item1_Col),
         dev_rad = bmm::wrap(deviation * pi / 180),
         across(ends_with("_Col"), ~ bmm::wrap((.x- Item1_Col)*pi/180) , .names = "{.col}_rad"),
         across(ends_with("_Pos"), ~ abs(bmm::wrap(2*pi*((.x - Item1_Pos)/13))),.names = "{.col}_rad")) %>% 
  mutate(
    LureIdx1 = case_when(SetSize >= 2 ~ 1,
                         TRUE ~ 0),
    LureIdx2 = case_when(SetSize >= 3 ~ 1,
                         TRUE ~ 0),
    LureIdx3 = case_when(SetSize >= 4 ~ 1,
                         TRUE ~ 0),
    LureIdx4 = case_when(SetSize >= 5 ~ 1,
                         TRUE ~ 0),
    LureIdx5 = case_when(SetSize >= 6 ~ 1,
                         TRUE ~ 0),
    LureIdx6 = case_when(SetSize >= 7 ~ 1,
                         TRUE ~ 0),
    LureIdx7 = case_when(SetSize >= 8 ~ 1, 
                         TRUE ~ 0),
    inv_SS = 1/(SetSize - 1),
    SetSize = as.factor(SetSize))


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
                        von_mises(link = "identity"),
                        von_mises(link = "identity"),
                        von_mises(link = "identity"),
                        order = "none")

# set up mixture model
IMM_mixModel_formula <- bf(dev_rad ~ 1,
                           # fix kappa over memory distributions
                           nlf(kappa1 ~ kappa),     # target distribution
                           nlf(kappa2 ~ kappa),     # non-target
                           nlf(kappa3 ~ kappa),     # non-target
                           nlf(kappa4 ~ kappa),     # non-target
                           nlf(kappa5 ~ kappa),     # non-target
                           nlf(kappa6 ~ kappa),     # non-target
                           nlf(kappa7 ~ kappa),     # non-target
                           nlf(kappa8 ~ kappa),     # non-target
                           # kappa for guessing distribution will be fixed using priors
                           kappa9 ~ 1,             # uniform  
                           # specify mixing distributions for distinct item categories
                           nlf(theta1 ~ exp(-s*Item1_Pos_rad)*c + a),   # p_mem
                           nlf(theta2 ~ LureIdx1*(exp(-s*Item2_Pos_rad)*c + a) + (1-LureIdx1)*(-100)),  # p_intrusion
                           nlf(theta3 ~ LureIdx2*(exp(-s*Item3_Pos_rad)*c + a) + (1-LureIdx2)*(-100)),  # p_intrusion
                           nlf(theta4 ~ LureIdx3*(exp(-s*Item4_Pos_rad)*c + a) + (1-LureIdx3)*(-100)),  # p_intrusion
                           nlf(theta5 ~ LureIdx4*(exp(-s*Item5_Pos_rad)*c + a) + (1-LureIdx4)*(-100)),  # p_intrusion
                           nlf(theta6 ~ LureIdx5*(exp(-s*Item6_Pos_rad)*c + a) + (1-LureIdx5)*(-100)),  # p_intrusion
                           nlf(theta7 ~ LureIdx6*(exp(-s*Item7_Pos_rad)*c + a) + (1-LureIdx6)*(-100)),  # p_intrusion
                           nlf(theta8 ~ LureIdx7*(exp(-s*Item8_Pos_rad)*c + a) + (1-LureIdx7)*(-100)),  # p_intrusion
                           # target & guessing distribution will be centered using priors
                           mu1 ~ 1, # fixed intercept constrained using priors
                           mu9 ~ 1, # fixed intercept constrained using priors
                           # center non-target distribution on data specified locations
                           nlf(mu2 ~ Item2_Col_rad),           # center non-target
                           nlf(mu3 ~ Item3_Col_rad),           # center non-target
                           nlf(mu4 ~ Item4_Col_rad),           # center non-target
                           nlf(mu5 ~ Item5_Col_rad),           # center non-target
                           nlf(mu6 ~ Item6_Col_rad),           # center non-target
                           nlf(mu7 ~ Item7_Col_rad),           # center non-target
                           nlf(mu8 ~ Item8_Col_rad),           # center non-target
                           nlf(s ~ exp(logS)),      
                           # now predict parameters of interest
                           kappa ~ 1 + SetSize + (1 || ID),  # fixed intercept & random slope: precision of memory distributions
                           logS ~ 1 + SetSize  + (1 || ID),   # fixed intercept & random slope: spatial gradient (on logarithmic scale)
                           c ~ 1 + SetSize + (1 || ID),      # fixed intercept & random slope: context activation
                           a ~ 1 + SetSize + (1 || ID),      # fixed intercept & random slope: general activation
                           # for brms to process this formula correctly, set non-linear to TRUE
                           nl = TRUE)

# check default priors
get_prior(IMM_mixModel_formula, df_OberauerLin2017_E1, IMM_mixModel)

# constrain priors to identify the model
IMM_priors <- 
  # first we center the target and guessing distribution to zero
  prior(constant(0), class = Intercept, dpar = "mu1") + 
  prior(constant(0), class = Intercept, dpar = "mu9") +
  # next, we set the guessing distribution to be uniform, kappa -> 0
  prior(constant(-100), class = Intercept, dpar = "kappa9") +
  # next, we set reasonable priors for the to be estimated distributions
  prior(normal(5.0, 0.8), class = b, nlpar = "kappa") +
  prior(normal(0, 1), class = b, nlpar = "logS", ) +
  prior(normal(2, 1), class = b, nlpar = "c") +
  prior(normal(0.5, 0.2), class = b, nlpar = "a")

# fir IMM using the brm function
fit_IMM_mixMod <- brm(formula = IMM_mixModel_formula, 
                      data = df_OberauerLin2017_E1,
                      family = IMM_mixModel, 
                      prior = IMM_priors, 
                      iter = 100 + 100, 
                      chains = nChains,
                      save_pars = save_pars(all = T))

# print out summary of results
summary(fit_IMM_mixMod)


cor(randFX$ID[,"Estimate","c_Intercept"],parms$c)
cor(randFX$ID[,"Estimate","a_Intercept"],parms$a)
cor(randFX$ID[,"Estimate","logS_Intercept"],log(parms$s))
cor(randFX$ID[,"Estimate","kappa_Intercept"],parms$kappa)
