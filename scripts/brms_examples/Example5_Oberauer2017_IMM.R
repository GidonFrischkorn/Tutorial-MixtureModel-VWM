#' This is the tutorial script for setting up the full Interference Measurement
#' model for visual working memory tasks that use continuous report recall procedures.
#' 
#' In this script, you will see:
#'  1) how the model is set up using the brms package, 
#'  2) how a simple version of the model is estimated, and 
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
warmup_samples <- 2000
postwarmup_samples <- 2000

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
                           theta9 ~ b,
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
                           kappa ~ 0 + SetSize + (0 + SetSize || ID),  # fixed intercept & random slope: precision of memory distributions
                           logS ~ 1,   # fixed intercept & random slope: spatial gradient (on logarithmic scale)
                           c ~ 1,      # fixed intercept & random slope: context activation
                           a ~ 0 + SetSize + (0 + SetSize || ID),      # fixed intercept & random slope: general activation
                           b ~ 0 + SetSize + (0 + SetSize || ID),
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
  prior(normal(2,2), class = b, nlpar = "c") +
  prior(constant(1), class = b, nlpar = "c") +
  prior(normal(1.5, 2), class = b, nlpar = "kappa") +
  prior(normal(0, 1), class = b, nlpar = "logS", ) +
  prior(normal(0,2), class = b, dpar = "theta9") +
  prior(normal(0.5, 1), class = b, nlpar = "a") +
  prior(constant(0), class = b, nlpar = "logS", coef = "SetSize1") +
  prior(constant(0), class = b, nlpar = "a", coef = "SetSize1")

if (!file.exists(here("output","fit_E5_OL2017_IMMfull.RData"))) {
  # fit IMM using the brm function
  fit_IMM_mixMod <- brm(formula = IMM_mixModel_formula, 
                        data = df_OberauerLin2017_E1,
                        family = IMM_mixModel, 
                        prior = IMM_priors,
                        
                        # save settings
                        sample_prior = TRUE,
                        save_pars = save_pars(all = TRUE),
                        
                        # add brms settings
                        warmup = warmup_samples,
                        iter = warmup_samples + postwarmup_samples, 
                        chains = nChains,
                        
                        # control commands for the sampler
                        control = list(adapt_delta = adapt_delta, 
                                       max_treedepth = max_treedepth))
  
  save(fit_IMM_mixMod,
       file = here("output","fit_E5_OL2017_IMMfull.RData"),
       compress = "xz")
} else {
  load(here("output","fit_E5_OL2017_IMMfull.RData"))
}

###############################################################################!
# 3) Model evaluation ----------------------------------------------------------
###############################################################################!

# plot the posterior predictive check to evaluate overall model fit
pp_check(fit_IMM_mixMod)

# print out summary of results
summary(fit_IMM_mixMod)

## 3.2) extract parameter estimates --------------------------------------------

# extract the fixed effects from the model
fixedEff <- fixef(fit_IMM_mixMod)

# determine the rows that contain the relevant parameter estimates
c_rows <- grepl("c_",rownames(fixedEff))
a_rows <- grepl("a_",rownames(fixedEff))
s_rows <- grepl("logS_",rownames(fixedEff))
kappa_rows <- grepl("kappa_",rownames(fixedEff))

# extract kappa estimates
c_fixedFX <- fixedEff[c_rows,]
a_fixedFX <- fixedEff[a_rows,]
s_fixedFX <- fixedEff[s_rows,]
kappa_fixedFX <- fixedEff[kappa_rows,]

# transform s & kappa from logarithmic to absolute scale
s_fixedFX <- exp(s_fixedFX)
kappa_fixedFX <- exp(kappa_fixedFX)

# print out parameter estimates
kappa_fixedFX
exp(c_fixedFX)
exp(a_fixedFX)
s_fixedFX

## 3.3) plot parameter estimates -----------------------------------------------
# extract posterior draws for fixed effects on kappa & theta
fixedFX_draws <- fit_IMM_mixMod %>%
  tidy_draws() %>%
  select(starts_with("b_"),.chain,.iteration,.draw) %>%
  pivot_longer(cols = starts_with("b_"),
               names_to = "modelPar",
               values_to = "postSample") %>%
  mutate(par = str_split_i(modelPar,"_",2),
         setsize = str_split_i(modelPar,"_",3),
         setsize = str_remove(setsize, "SetSize")) %>%
  select(-modelPar) %>%
  filter(par %in% c("c","a","logS","kappa")) %>%
  mutate(postSample_abs = case_when(par %in% c("logS","kappa") ~ exp(postSample),
                                    TRUE ~ postSample))

# plot kappa results
kappa_plot <- ggplot(data = fixedFX_draws %>% filter(par == "kappa"),
                     aes(x = setsize, y = postSample_abs)) +
  coord_cartesian(ylim = c(0,30)) +
  geom_half_violin(position = position_nudge(x = .1, y = 0), side = "r", fill = "darkgrey", color = NA,
                   adjust = 1, trim = TRUE, alpha = 0.9, show.legend = FALSE, scale = "width") +
  stat_summary(geom = "pointrange", fun.data = mode_hdi,
               size = 0.3, linewidth = 0.8,
               position = position_dodge(0.1)) +
  # geom_point(data = results_LS_2018 %>% filter(param == "contSD"),
  #            aes(y = mean, x = RI, color = as.factor(nCues)),
  #            shape = "diamond", size = 2.5,
  #            position = position_nudge(x = -.1, y = 0)) +
  scale_fill_grey(start = 0, end = .8) +
  scale_color_grey(start = 0, end = .8) +
  labs(x = "Set Size", y = "Memory imprecision (SD)", fill = "No. of Cues", color = "No. of Cues",
       title = "B") +
  guides(color = "none") +
  clean_plot()
kappa_plot

# plot pMem results
c_plot <- ggplot(data = fixedFX_draws %>% filter(par == "c"),
                     aes(x = setsize, y = postSample_abs)) +
  coord_cartesian(ylim = c(0,8)) +
  geom_half_violin(position = position_nudge(x = .1, y = 0), side = "r", fill = "darkgrey", color = NA,
                   adjust = 1, trim = TRUE, alpha = 0.9, show.legend = FALSE, scale = "width") +
  stat_summary(geom = "pointrange", fun.data = mode_hdi,
               size = 0.3, linewidth = 0.8,
               position = position_dodge(0.1)) +
  # geom_point(data = results_LS_2018 %>% filter(param == "contSD"),
  #            aes(y = mean, x = RI, color = as.factor(nCues)),
  #            shape = "diamond", size = 2.5,
  #            position = position_nudge(x = -.1, y = 0)) +
  scale_fill_grey(start = 0, end = .8) +
  scale_color_grey(start = 0, end = .8) +
  labs(x = "Set Size", y = "Context Activation (c)",
       title = "B") +
  guides(color = "none") +
  clean_plot()
c_plot

# plot pMem results
a_plot <- ggplot(data = fixedFX_draws %>% filter(par == "a", setsize != "1"),
                 aes(x = setsize, y = postSample_abs)) +
  coord_cartesian(ylim = c(-4,1)) +
  geom_half_violin(position = position_nudge(x = .1, y = 0), side = "r", fill = "darkgrey", color = NA,
                   adjust = 1, trim = TRUE, alpha = 0.9, show.legend = FALSE, scale = "width") +
  stat_summary(geom = "pointrange", fun.data = mode_hdi,
               size = 0.3, linewidth = 0.8,
               position = position_dodge(0.1)) +
  # geom_point(data = results_LS_2018 %>% filter(param == "contSD"),
  #            aes(y = mean, x = RI, color = as.factor(nCues)),
  #            shape = "diamond", size = 2.5,
  #            position = position_nudge(x = -.1, y = 0)) +
  scale_fill_grey(start = 0, end = .8) +
  scale_color_grey(start = 0, end = .8) +
  labs(x = "Set Size", y = "General Activation (a)",
       title = "B") +
  guides(color = "none") +
  clean_plot()
a_plot

# plot pMem results
s_plot <- ggplot(data = fixedFX_draws %>% filter(par == "logS", setsize != "1"),
                 aes(x = setsize, y = postSample_abs)) +
  coord_cartesian(ylim = c(-10,150)) +
  geom_half_violin(position = position_nudge(x = .1, y = 0), side = "r", fill = "darkgrey", color = NA,
                   adjust = 1.5, trim = TRUE, alpha = 0.9, show.legend = FALSE, scale = "width") +
  stat_summary(geom = "pointrange", fun.data = mode_hdi,
               size = 0.3, linewidth = 0.8,
               position = position_dodge(0.1)) +
  # geom_point(data = results_LS_2018 %>% filter(param == "contSD"),
  #            aes(y = mean, x = RI, color = as.factor(nCues)),
  #            shape = "diamond", size = 2.5,
  #            position = position_nudge(x = -.1, y = 0)) +
  scale_fill_grey(start = 0, end = .8) +
  scale_color_grey(start = 0, end = .8) +
  labs(x = "Set Size", y = "Spatial Specificify (s)",
       title = "B") +
  guides(color = "none") +
  clean_plot()
s_plot

# patch plots together
joint_plot <- kappa_plot+ a_plot + c_plot + s_plot + 
  plot_layout(ncol = 2)

# show joint plot
joint_plot

# save plots with high resolution
ggsave(
  filename = here("figures","plot_kappaEst_OL2017.jpeg"),
  plot = kappa_plot, width = 6, height = 6
)

ggsave(
  filename = here("figures","plot_cEst_OL2017.jpeg"),
  plot = c_plot, width = 6, height = 6
)

ggsave(
  filename = here("figures","plot_aEst_OL2017.jpeg"),
  plot = a_plot, width = 6, height = 6
)

ggsave(
  filename = here("figures","plot_sEst_OL2017.jpeg"),
  plot = s_plot, width = 6, height = 6
)

ggsave(
  filename = here("figures","plot_jointRes_LS2018.jpeg"),
  plot = joint_plot, width = 6*2, height = 6
)
