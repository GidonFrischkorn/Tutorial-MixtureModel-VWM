#' This is the tutorial script for setting up the Interference Measurement Model (IMMbsc)
#' for visual working memory tasks that use continuous report recall procedures.
#' 
#' In this script, you will see:
#'  1) how the model is set up using the bmm package, 
#'  2) how a simple version of the model is estimates, and 
#'  3) how the model can be evaluated and results extracted and plotted.

# 0) R Setup: Packages & Data --------------------------------------------------

# load required packages
pacman::p_load(here, bmm, brms, tidyverse, tidybayes, patchwork, gghalves)

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
data <- oberauer_lin_2017

###############################################################################!
# Fit IMMbsc ----------------------------------------------------------------
###############################################################################!

imm_bsc_model <- imm(resp_error = "dev_rad",
                     nt_features = "col_nt", 
                     nt_distances = "dist_nt", regex = TRUE,
                     set_size = "set_size",
                     version = "bsc")

# set up mixture model
imm_bsc_formula <- bmf(# fixed intercept & random slope: precision of memory distributions
         kappa ~ 0 + set_size + (0 + set_size || ID),
         # fixed intercept & random slope: context activation
         c ~ 0 + set_size + (0 + set_size || ID),
         # fixed intercept & random slope: general activation
         s ~ 0 + set_size + (0 + set_size || ID))

filename_IMMbsc <- here("output","fit_E5_OL2017_IMMbsc")

imm_bsc_fit <- bmm(
  model = imm_bsc_model,
  formula = imm_bsc_formula, 
  data = data,
  
  # save settings
  sample_prior = TRUE,
  save_pars = save_pars(all = TRUE),
  
  # add brms settings
  warmup = warmup_samples,
  iter = warmup_samples + postwarmup_samples, 
  chains = nChains,
  cores = parallel::detectCores(),
  
  # control commands for the sampler
  control = list(adapt_delta = adapt_delta, 
                 max_treedepth = max_treedepth),
  
  # save results to file
  file = filename_IMMbsc
)

imm_bsc_fit$formula

###############################################################################!
# Model evaluation ----------------------------------------------------------
###############################################################################!

# plot the posterior predictive check to evaluate overall model fit
pp_check(imm_bsc_fit, group = "set_size", type = "dens_overlay_grouped")

# print out summary of results
summary(imm_bsc_fit)

## extract parameter estimates --------------------------------------------

# extract the fixed effects from the model
fixedEff <- fixef(imm_bsc_fit)

# determine the rows that contain the relevant parameter estimates
c_rows <- grepl("c_",rownames(fixedEff))
s_rows <- grepl("s_",rownames(fixedEff))
kappa_rows <- grepl("kappa_",rownames(fixedEff))

# extract kappa estimates
c_fixedFX <- fixedEff[c_rows,]
s_fixedFX <- fixedEff[s_rows,]
kappa_fixedFX <- fixedEff[kappa_rows,]

# transform s & kappa from logarithmic to absolute scale
kappa_fixedFX <- exp(kappa_fixedFX)

# print out parameter estimates
kappa_fixedFX
exp(c_fixedFX)
exp(s_fixedFX)

## plot parameter estimates -----------------------------------------------

# extract posterior draws for fixed effects on kappa & theta
fixedFX_draws <- imm_bsc_fit %>%
  tidy_draws() %>%
  select(starts_with("b_"),.chain,.iteration,.draw) %>%
  pivot_longer(cols = starts_with("b_"),
               names_to = "modelPar",
               values_to = "postSample") %>%
  mutate(par = str_split_i(modelPar,"_",2),
         setsize = str_split_i(modelPar,"_",4),
         setsize = str_remove(setsize, "size")) %>%
  select(-modelPar) %>%
  filter(par %in% c("c","a","s","kappa")) %>%
  mutate(postSample_abs = case_when(par %in% c("c","s","kappa") ~ exp(postSample),
                                    TRUE ~ postSample))

# plot kappa results
plot_kappa_IMMbsc <- ggplot(data = fixedFX_draws %>% filter(par == "kappa"),
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
plot_c_IMMbsc <- ggplot(data = fixedFX_draws %>% filter(par == "c"),
                        aes(x = setsize, y = postSample_abs)) +
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
  scale_y_log10() +
  labs(x = "Set Size", y = "Context Activation (c)",
       title = "A") +
  guides(color = "none") +
  clean_plot()

# plot pMem results
plot_s_IMMbsc <- ggplot(data = fixedFX_draws %>% filter(par == "s", setsize != "1"),
                        aes(x = setsize, y = (postSample_abs))) +
  coord_cartesian(ylim = c(0,40)) +
  geom_half_violin(position = position_nudge(x = .1, y = 0), 
                   side = "r", fill = "darkgrey", color = NA,
                   adjust = 1, trim = TRUE, alpha = 0.9, 
                   show.legend = FALSE, scale = "width") +
  stat_summary(geom = "pointrange", fun.data = mode_qi,
               size = 0.3, linewidth = 0.8,
               position = position_dodge(0.1)) +
  # geom_point(data = results_LS_2018 %>% filter(param == "contSD"),
  #            aes(y = mean, x = RI, color = as.factor(nCues)),
  #            shape = "diamond", size = 2.5,
  #            position = position_nudge(x = -.1, y = 0)) +
  scale_fill_grey(start = 0, end = .8) +
  scale_color_grey(start = 0, end = .8) +
  labs(x = "Set Size", y = "Generalization Gradient (s)",
       title = "B") +
  guides(color = "none") +
  clean_plot()


# patch plots together
joint_plot <-   plot_c_IMMbsc + plot_s_IMMbsc + plot_kappa_IMMbsc +
  plot_layout(ncol = 3)

# show joint plot
joint_plot

# save plots with high resolution
ggsave(
  filename = here("figures","plotAll_OL2017_IMMbsc.jpeg"),
  plot = joint_plot, width = 4*3, height = 4
)
