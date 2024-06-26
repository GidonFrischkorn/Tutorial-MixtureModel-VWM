#' This is the tutorial script for setting up the Interference Measurement Model (IMMabc)
#' for visual working memory tasks that use continuous report recall procedures for the 
#' data from Oberauer & Lin (2017)
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
postwarmup_samples<- 2000

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
# Fit IMMabc ----------------------------------------------------------------
###############################################################################!

imm_abc_model <- imm(resp_error = "dev_rad",
                     nt_features = "col_nt", regex = TRUE,
                     set_size = "set_size",
                     version = "abc")

## Estimate pars for IMMabc ------------------------------------------------
# set up mixture model
imm_abc_formula <- bmf(
  # fixed intercept & random slope: precision of memory distributions
  kappa ~ 0 + set_size + (0 + set_size || ID),
  # fixed intercept & random slope: context activation
  c ~ 0 + set_size + (0 + set_size || ID),
  # fixed intercept & random slope: general activation
  a ~ 0 + set_size + (0 + set_size || ID)
)

filename_IMMabc = here("output","fit_E5_OL2017_IMMabc")

# fit IMM using the bmm function
imm_abc_fit <- bmm(
  model = imm_abc_model,
  formula = imm_abc_formula, 
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
  file = filename_IMMabc
)

# plot the posterior predictive check to evaluate overall model fit
pp_check(imm_abc_fit, group = "set_size", type = "dens_overlay_grouped")

# print out summary of results
summary(imm_abc_fit)

###############################################################################!
# Model evaluation ----------------------------------------------------------
###############################################################################!

## extract parameter estimates --------------------------------------------

# extract the fixed effects from the model
fixedEff <- fixef(imm_abc_fit)

# determine the rows that contain the relevant parameter estimates
c_rows <- grepl("c_",rownames(fixedEff))
a_rows <- grepl("a_",rownames(fixedEff)) & !grepl("kappa_",rownames(fixedEff))
kappa_rows <- grepl("kappa_",rownames(fixedEff))

# extract kappa estimates
c_fixedFX <- fixedEff[c_rows,]
a_fixedFX <- fixedEff[a_rows,]
kappa_fixedFX <- fixedEff[kappa_rows,]

# transform s & kappa from logarithmic to absolute scale
kappa_fixedFX <- exp(kappa_fixedFX)

# print out parameter estimates
kappa_fixedFX
exp(c_fixedFX)
exp(a_fixedFX)

## plot parameter estimates -----------------------------------------------

# extract posterior draws for fixed effects on kappa & theta
fixedFX_draws <- imm_abc_fit %>%
  tidy_draws() %>%
  select(starts_with("b_"),.chain,.iteration,.draw) %>%
  pivot_longer(cols = starts_with("b_"),
               names_to = "modelPar",
               values_to = "postSample") %>%
  mutate(par = str_split_i(modelPar,"_",2),
         setsize = str_split_i(modelPar,"_",4),
         setsize = str_remove(setsize, "size")) %>%
  select(-modelPar) %>%
  filter(par %in% c("c","a","kappa")) %>%
  mutate(postSample_abs = case_when(par %in% c("c","a","kappa") ~ exp(postSample),
                                    TRUE ~ postSample))

df_samples_probConversion <- fixedFX_draws %>% 
  select(-postSample) %>% 
  pivot_wider(values_from = postSample_abs,
              names_from = par) %>%
  mutate(setsize = as.numeric(setsize),
         Pmem = (c + a)/((c + a) + (setsize-1)*(a) + exp(0)),
         Pnt = case_when(setsize == 1 ~ 0,
                         TRUE ~ ((setsize-1)*(a))/((c + a) + (setsize-1)*(a) + exp(0))),
         Pg = 1 - Pmem - Pnt,
         Ptotal = Pmem + Pnt + Pg)

# plot derived three-parameter mixture parameters
plot_Pmem <- ggplot(data = df_samples_probConversion,
       aes(y = Pmem, x = as.factor(setsize))) +
  coord_cartesian(ylim = c(0,1)) +
  geom_half_violin(position = position_nudge(x = .1, y = 0), side = "r", fill = "darkgrey", color = NA,
                                                                adjust = 1, trim = TRUE, alpha = 0.9, show.legend = FALSE, scale = "width") +
  stat_summary(geom = "pointrange", fun.data = mode_hdi,
               size = 0.3, linewidth = 0.8,
               position = position_dodge(0.1)) +
  scale_fill_grey(start = 0, end = .8) +
  scale_color_grey(start = 0, end = .8) +
  labs(x = "Set Size", y = "Pmem", fill = "No. of Cues", color = "No. of Cues",
       title = "A") +
  guides(color = "none") +
  clean_plot()

plot_Pnt <- ggplot(data = df_samples_probConversion,
       aes(y = Pnt, x = as.factor(setsize))) +
  coord_cartesian(ylim = c(0,1)) +
  geom_half_violin(position = position_nudge(x = .1, y = 0), side = "r", fill = "darkgrey", color = NA,
                   adjust = 1, trim = TRUE, alpha = 0.9, show.legend = FALSE, scale = "width") +
  stat_summary(geom = "pointrange", fun.data = mode_hdi,
               size = 0.3, linewidth = 0.8,
               position = position_dodge(0.1)) +
  scale_fill_grey(start = 0, end = .8) +
  scale_color_grey(start = 0, end = .8) +
  labs(x = "Set Size", y = "Pnt", fill = "No. of Cues", color = "No. of Cues",
       title = "B") +
  guides(color = "none") +
  clean_plot()

plot_Pg <- ggplot(data = df_samples_probConversion,
       aes(y = Pg, x = as.factor(setsize))) +
  coord_cartesian(ylim = c(0,1)) +
  geom_half_violin(position = position_nudge(x = .1, y = 0), side = "r", fill = "darkgrey", color = NA,
                   adjust = 1, trim = TRUE, alpha = 0.9, show.legend = FALSE, scale = "width") +
  stat_summary(geom = "pointrange", fun.data = mode_hdi,
               size = 0.3, linewidth = 0.8,
               position = position_dodge(0.1)) +
  scale_fill_grey(start = 0, end = .8) +
  scale_color_grey(start = 0, end = .8) +
  labs(x = "Set Size", y = "Pg", fill = "No. of Cues", color = "No. of Cues",
       title = "C") +
  guides(color = "none") +
  clean_plot()

plot_kappa <- ggplot(data = df_samples_probConversion,
       aes(y = kappa, x = as.factor(setsize))) +
  coord_cartesian(ylim = c(0,30)) +
  geom_half_violin(position = position_nudge(x = .1, y = 0), side = "r", fill = "darkgrey", color = NA,
                   adjust = 1, trim = TRUE, alpha = 0.9, show.legend = FALSE, scale = "width") +
  stat_summary(geom = "pointrange", fun.data = mode_hdi,
               size = 0.3, linewidth = 0.8,
               position = position_dodge(0.1)) +
  scale_fill_grey(start = 0, end = .8) +
  scale_color_grey(start = 0, end = .8) +
  labs(x = "Set Size", y = "Memory precision (kappa)", fill = "No. of Cues", color = "No. of Cues",
       title = "D") +
  guides(color = "none") +
  clean_plot()

joint_plot <-   plot_Pmem + plot_Pnt + plot_Pg + plot_kappa +
  plot_layout(ncol = 2, nrow = 2)
joint_plot

ggsave(
  filename = here("figures","plotAll_OL2017_IMMabc_3pPars.jpeg"),
  plot = joint_plot, width = 6, height = 6
)

# plot kappa results
plot_kappa_IMMabc <- ggplot(data = fixedFX_draws %>% filter(par == "kappa"),
                            aes(x = setsize, y = postSample_abs)) +
  coord_cartesian(ylim = c(0,30)) +
  geom_half_violin(position = position_nudge(x = .1, y = 0), side = "r", fill = "darkgrey", color = NA,
                   adjust = 1, trim = TRUE, alpha = 0.9, show.legend = FALSE, scale = "width") +
  stat_summary(geom = "pointrange", fun.data = mode_hdi,
               size = 0.3, linewidth = 0.8,
               position = position_dodge(0.1)) +
  scale_fill_grey(start = 0, end = .8) +
  scale_color_grey(start = 0, end = .8) +
  labs(x = "Set Size", y = "Memory precision (kappa)", fill = "No. of Cues", color = "No. of Cues",
       title = "C") +
  guides(color = "none") +
  clean_plot()

# plot pMem results
plot_c_IMMabc <- ggplot(data = fixedFX_draws %>% filter(par == "c"),
                        aes(x = setsize, y = postSample_abs)) +
  #coord_cartesian(ylim = c(0,10)) +
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
plot_a_IMMabc <- ggplot(data = fixedFX_draws %>% filter(par == "a", setsize != "1"),
                        aes(x = setsize, y = postSample_abs)) +
  coord_cartesian(ylim = c(0,2.5)) +
  geom_hline(yintercept = 1, color ="firebrick",
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
