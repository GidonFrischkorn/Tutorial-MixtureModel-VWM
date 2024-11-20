#' This is the tutorial script for setting up the Three-parameter mixture model (mixture3p)
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

#############################################################################!
# 1) Fit 3par Mixture model - using bmm functions                         ####
#############################################################################!

## Estimate pars for 3par Mixture Model ------------------------------------
model_3pMM <- mixture3p(resp_error = "dev_rad",
                        nt_features = paste0("col_nt",1:7),
                        set_size = "set_size")

# formula
formula_3pMM <- bmf(kappa ~ 0 + set_size + (0 + set_size || ID),
                    thetat ~ 0 + set_size + (0 + set_size || ID),
                    thetant ~ 0 + set_size + (0 + set_size || ID))

default_prior(ff, model = model_3p, data = data)

# if the model has been already estimated, load the results, otherwise estimate it
filename = here("output","fit_E5_OL2017_3pMM")

fit_3pMM <- bmm(
  model = model_3p,
  formula = formula_3pMM, 
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
  
  # save results in file
  file = filename
)

# plot the posterior predictive check to evaluate overall model fit
pp_check(fit_3pMM, group = "set_size", type = "dens_overlay_grouped")

# print out summary of results
summary(fit_3pMM)


###############################################################################!
# 3) Plot 3par Mixture Model results -----------------------------------------------
###############################################################################!

## Plot 3par Mixture Model results -----------------------------------------
# extract posterior estimates from the model
draws_3pMM <- fit_3pMM %>%
  tidy_draws() %>%
  select(starts_with("b_"),.chain,.iteration,.draw) %>%
  pivot_longer(cols = starts_with("b_"),
               names_to = "modelPar",
               values_to = "postSample") %>%
  mutate(par = str_split_i(modelPar,"_",2),
         setsize = str_split_i(modelPar,"_",4),
         setsize = str_remove(setsize, "size")) %>%
  select(-modelPar) %>%
  filter(par %in% c("thetat","thetant","kappa")) %>%
  mutate(postSample_abs = case_when(par %in% c("kappa") ~ exp(postSample),
                                    TRUE ~ postSample))

plot_kappa_3pMM <- ggplot(data = draws_3pMM %>% filter(par == "kappa"),
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
       title = "D") +
  guides(color = "none") +
  clean_plot()

# convert mixture weights into recall probabilities
thetaSamples_3pMM <- draws_3pMM %>% 
  dplyr::filter(par != "kappa") %>% 
  dplyr::select(.chain,.iteration,.draw,setsize,par,postSample) %>% 
  pivot_wider(names_from = par,
              values_from = postSample) %>% 
  mutate(p_Mem = exp(thetat)/(exp(thetat) + exp(thetant) + exp(0)),
         p_Swap = exp(thetant)/(exp(thetat) + exp(thetant) + exp(0)),
         p_Guess = exp(0)/(exp(thetat) + exp(thetant) + exp(0)),
  ) %>% 
  dplyr::select(.chain,.iteration,.draw,setsize,p_Mem,p_Swap,p_Guess) %>% 
  pivot_longer(cols = starts_with("p_"),
               names_to = "coef",
               values_to = "sample") 

plot_pMem_3pMM <- ggplot(data = thetaSamples_3pMM %>% filter(coef == "p_Mem"),
                         aes(x = setsize, y = sample)) +
  coord_cartesian(ylim = c(0,1)) +
  geom_half_violin(position = position_nudge(x = .1, y = 0), side = "r", fill = "darkgrey", color = NA,
                   adjust = 1, trim = TRUE, alpha = 0.9, show.legend = FALSE, scale = "width") +
  stat_summary(geom = "pointrange", fun.data = mode_hdi,
               size = 0.3, linewidth = 0.8,
               position = position_dodge(0.1)) +
  scale_fill_grey(start = 0, end = .8) +
  scale_color_grey(start = 0, end = .8) +
  labs(x = "Set Size", y = expression(P[mem]), fill = "No. of Cues", color = "No. of Cues",
       title = "A") +
  guides(color = "none") +
  clean_plot()

plot_pSwap_3pMM <- ggplot(data = thetaSamples_3pMM %>% filter(coef == "p_Swap"),
                          aes(x = setsize, y = sample)) +
  coord_cartesian(ylim = c(0,1)) +
  geom_half_violin(position = position_nudge(x = .1, y = 0), side = "r", fill = "darkgrey", color = NA,
                   adjust = 1, trim = TRUE, alpha = 0.9, show.legend = FALSE, scale = "width") +
  stat_summary(geom = "pointrange", fun.data = mode_hdi,
               size = 0.3, linewidth = 0.8,
               position = position_dodge(0.1)) +
  scale_fill_grey(start = 0, end = .8) +
  scale_color_grey(start = 0, end = .8) +
  labs(x = "Set Size", y = expression(P[swap]), fill = "No. of Cues", color = "No. of Cues",
       title = "B") +
  guides(color = "none") +
  clean_plot()

plot_pGuess_3pMM <- ggplot(data = thetaSamples_3pMM %>% filter(coef == "p_Guess"),
                           aes(x = setsize, y = sample)) +
  coord_cartesian(ylim = c(0,1)) +
  geom_half_violin(position = position_nudge(x = .1, y = 0), side = "r", fill = "darkgrey", color = NA,
                   adjust = 1, trim = TRUE, alpha = 0.9, show.legend = FALSE, scale = "width") +
  stat_summary(geom = "pointrange", fun.data = mode_hdi,
               size = 0.3, linewidth = 0.8,
               position = position_dodge(0.1)) +
  scale_fill_grey(start = 0, end = .8) +
  scale_color_grey(start = 0, end = .8) +
  labs(x = "Set Size", y = expression(P[guess]), fill = "No. of Cues", color = "No. of Cues",
       title = "C") +
  guides(color = "none") +
  clean_plot()

joint_plot_3pMM <- plot_pMem_3pMM + plot_pSwap_3pMM + plot_pGuess_3pMM + plot_kappa_3pMM +
  plot_layout(ncol = 2, nrow = 2)
joint_plot_3pMM

ggsave(
  filename = here("figures","plotAll_OL2017_3pMM.jpeg"),
  plot = joint_plot_3pMM, width = 6, height = 6
)
