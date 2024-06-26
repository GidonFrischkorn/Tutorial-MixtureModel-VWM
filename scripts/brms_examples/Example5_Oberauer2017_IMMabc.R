#' This is the tutorial script for setting up the Interference Measurement Model
#' in the abc version, assuming that swaps occur independent of spatial proximity
#' to the target for visual working memory tasks that use continuous report recall procedures.
#' 
#' In this script, you will see:
#'  1) how the model is set up using the brms package, 
#'  2) how a simple version of the model is estimated, and 
#'  3) how the model can be evaluated and results extracted and plotted.

# 0) R Setup: Packages & Data --------------------------------------------------

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

# compute relevant variables for estimating the 3-parameter mixture modelpo
df_OberauerLin2017_E1 <- df_OberauerLin2017_E1 %>% 
  mutate(deviation = (Response - Item1_Col),
         devRad = bmm::wrap(deviation * pi / 180),
         across(ends_with("_Col"), ~ bmm::wrap((.x- Item1_Col)*pi/180) , .names = "{.col}_rad"),
         across(ends_with("_Pos"), ~ abs(bmm::wrap(2*pi*((.x - Item1_Pos)/13))),.names = "{.col}_rad"),
         SetSize = as.factor(SetSize))

#############################################################################!
# 1) Fit 3par Mixture model - using bmm functions                         ####
#############################################################################!

# set up the model object
mixture3p_model <- mixture3p(resp_error = "devRad",
                             nt_features = paste0("Item",2:8,"_Col_rad"),
                             set_size = "SetSize")

# formula
ff <- bmf(kappa ~ 0 + SetSize + (0 + SetSize || ID),
          thetat ~ 0 + SetSize + (0 + SetSize || ID),
          thetant ~ 0 + SetSize + (0 + SetSize || ID))

# if the model has been already estimated, load the results, otherwise estimate it
filename <- "output/fit_E5_OL2017_3pMM"

fit_3pMM <- bmm::bmm(
  formula = ff, 
  data = df_OberauerLin2017_E1, 
  model = mixture3p_model,
  
  # save settings
  sample_prior = TRUE,
  save_pars = save_pars(all = TRUE),
  
  # add brms settings
  warmup = warmup_samples,
  iter = warmup_samples + postwarmup_samples, 
  chains = nChains,
  cores = nChains,
  
  # control commands for the sampler
  control = list(adapt_delta = adapt_delta, 
                 max_treedepth = max_treedepth),
  
  file = filename
)

# plot the posterior predictive check to evaluate overall model fit
pp_check(fit_3pMM)

# print out summary of results
summary(fit_3pMM)


###############################################################################!
# 2) Fit IMMabc ----------------------------------------------------------------
###############################################################################!

df_OberauerLin2017_E1 <- df_OberauerLin2017_E1 %>% 
  mutate(
    LureIdx1 = case_when(as.numeric(SetSize) >= 2 ~ 1,
                         TRUE ~ 0),
    LureIdx2 = case_when(as.numeric(SetSize)  >= 3 ~ 1,
                         TRUE ~ 0),
    LureIdx3 = case_when(as.numeric(SetSize)  >= 4 ~ 1,
                         TRUE ~ 0),
    LureIdx4 = case_when(as.numeric(SetSize)  >= 5 ~ 1,
                         TRUE ~ 0),
    LureIdx5 = case_when(as.numeric(SetSize)  >= 6 ~ 1,
                         TRUE ~ 0),
    LureIdx6 = case_when(as.numeric(SetSize)  >= 7 ~ 1,
                         TRUE ~ 0),
    LureIdx7 = case_when(as.numeric(SetSize)  >= 8 ~ 1, 
                         TRUE ~ 0),
    inv_SS = 1/(as.numeric(SetSize)  - 1))

# create mixture of von Mises distributions
IMMabc_mixFamily <- mixture(von_mises(link = "identity"),
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
IMMabc_mixModel_formula <- bf(devRad ~ 1,
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
                              nlf(theta1 ~ log(exp(c) + exp(a))),   # p_mem
                              nlf(theta2 ~ LureIdx1*(a) + (1-LureIdx1)*(-100)),  # p_intrusion
                              nlf(theta3 ~ LureIdx2*(a) + (1-LureIdx2)*(-100)),  # p_intrusion
                              nlf(theta4 ~ LureIdx3*(a) + (1-LureIdx3)*(-100)),  # p_intrusion
                              nlf(theta5 ~ LureIdx4*(a) + (1-LureIdx4)*(-100)),  # p_intrusion
                              nlf(theta6 ~ LureIdx5*(a) + (1-LureIdx5)*(-100)),  # p_intrusion
                              nlf(theta7 ~ LureIdx6*(a) + (1-LureIdx6)*(-100)),  # p_intrusion
                              nlf(theta8 ~ LureIdx7*(a) + (1-LureIdx7)*(-100)),  # p_intrusion
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
                              # now predict parameters of interest
                              kappa ~ 0 + SetSize + (0 + SetSize || ID),  # fixed intercept & random slope: precision of memory distributions
                              c ~ 0 + SetSize + (0 + SetSize || ID),      # fixed intercept & random slope: context activation
                              a ~ 0 + SetSize + (0 + SetSize || ID),      # fixed intercept & random slope: general activation
                              # for brms to process this formula correctly, set non-linear to TRUE
                              nl = TRUE)

# check default priors
get_prior(IMMabc_mixModel_formula, df_OberauerLin2017_E1, IMMabc_mixFamily)

# constrain priors to identify the model
IMMabc_priors <- 
  # first we center the target and guessing distribution to zero
  prior(constant(0), class = Intercept, dpar = "mu1") + 
  prior(constant(0), class = Intercept, dpar = "mu9") +
  # next, we set the guessing distribution to be uniform, kappa -> 0
  prior(constant(-100), class = Intercept, dpar = "kappa9") +
  # next, we set priors for the to be estimated distributions
  prior(normal(0,1), class = b, nlpar = "c") +
  prior(normal(0, 2), class = b, nlpar = "kappa") +
  prior(normal(0, 1), class = b, nlpar = "a") +
  prior(constant(0), class = b, nlpar = "a", coef = "SetSize1")

filename <- "output/fit_E5_OL2017_IMMabc_brms"

# fit IMM using the brm function
fit_IMMabc_mixMod <- brm(formula = IMMabc_mixModel_formula, 
                         data = df_OberauerLin2017_E1,
                         family = IMMabc_mixFamily, 
                         prior = IMMabc_priors,
                         
                         # save settings
                         sample_prior = TRUE,
                         save_pars = save_pars(all = TRUE),
                         
                         # add brms settings
                         warmup = warmup_samples,
                         iter = warmup_samples + postwarmup_samples, 
                         chains = nChains,
                         
                         # control commands for the sampler
                         control = list(adapt_delta = adapt_delta, 
                                        max_treedepth = max_treedepth),
                         
                         file = filename)


###############################################################################!
# 3) Model evaluation ----------------------------------------------------------
###############################################################################!

# plot the posterior predictive check to evaluate overall model fit
pp_check(fit_IMMabc_mixMod)

# print out summary of results
summary(fit_IMMabc_mixMod)

## 3.2) extract parameter estimates --------------------------------------------

# extract the fixed effects from the model
fixedEff <- fixef(fit_IMMabc_mixMod)

# determine the rows that contain the relevant parameter estimates
c_rows <- grepl("c_",rownames(fixedEff))
a_rows <- grepl("a_",rownames(fixedEff)) & !grepl("kappa_",rownames(fixedEff))
s_rows <- grepl("logS_",rownames(fixedEff))
kappa_rows <- grepl("kappa_",rownames(fixedEff))

# extract kappa estimates
c_fixedFX <- fixedEff[c_rows,]
a_fixedFX <- fixedEff[a_rows,]
s_fixedFX <- fixedEff[s_rows,]
kappa_fixedFX <- fixedEff[kappa_rows,]

# transform s & kappa from logarithmic to absolute scale
kappa_fixedFX <- exp(kappa_fixedFX)

# print out parameter estimates
kappa_fixedFX
exp(c_fixedFX)
exp(a_fixedFX)

## 3.3) plot parameter estimates -----------------------------------------------
results_OL_2017 <- read.table(here("data","LS2018_2P_hierarchicalfit.txt"),
                              header = T, sep = ",") %>%
  filter(param != "catActive") %>%
  mutate(RI = retention,
         nCues = case_when(cueCond == "NoCue" ~ 0,
                           cueCond == "RetroCue" & RI == "short" ~ 1,
                           cueCond == "RetroCue" & RI == "long" ~ 2),
         ageGroup = case_when(BP_Group == "Old" ~ "Old",
                              TRUE ~ "Young"))

# extract posterior draws for fixed effects on kappa & theta
fixedFX_draws <- fit_IMMabc_mixMod %>%
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
plot_kappa_IMMabc <- ggplot(data = fixedFX_draws %>% filter(par == "kappa"),
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
  labs(x = "Set Size", y = "Memory precision (kappa)", fill = "No. of Cues", color = "No. of Cues",
       title = "C") +
  guides(color = "none") +
  clean_plot()

# plot pMem results
plot_c_IMMabc <- ggplot(data = fixedFX_draws %>% filter(par == "c"),
                        aes(x = setsize, y = exp(postSample_abs))) +
  coord_cartesian(ylim = c(0,100)) +
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
       title = "A") +
  guides(color = "none") +
  clean_plot()

# plot pMem results
plot_a_IMMabc <- ggplot(data = fixedFX_draws %>% filter(par == "a", setsize != "1"),
                        aes(x = setsize, y = exp(postSample_abs))) +
  coord_cartesian(ylim = c(0,2)) +
  geom_hline(yintercept = exp(0), color ="firebrick", 
             linetype = "dotted", linewidth = 1) +
  geom_half_violin(position = position_nudge(x = .1, y = 0), 
                   side = "r", fill = "darkgrey", color = NA,
                   adjust = 1, trim = TRUE, alpha = 0.9, 
                   show.legend = FALSE, scale = "width") +
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


# patch plots together
joint_plot <-   plot_c_IMMabc + plot_a_IMMabc + plot_kappa_IMMabc +
  plot_layout(ncol = 3)

# show joint plot
joint_plot

# save plots with high resolution
ggsave(
  filename = here("figures","plotAll_OL2017_IMMabc.jpeg"),
  plot = joint_plot, width = 4*3, height = 4
)
