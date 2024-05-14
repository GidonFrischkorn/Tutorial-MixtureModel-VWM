# Reproducing results from Loaiza & Souza (2018)
# "Is Refreshing in Working Memory Impaired in Older Age? Evidence from the Retro-Cue Paradigm"

# start fresh
rm(list = ls()) # clean up work space
graphics.off()  # switch off graphics device

# load required packages
pacman::p_load(here, brms, tidyverse, tidybayes, patchwork, bayestestR, gghalves)
pacman::p_load_gh("venpopov/bmm")

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
nChains <- 6

#' if the number of user defined chains is larger than the number of cores 
#' on the system than estimate as many chains as there are cores on the system
if (nChains >  parallel::detectCores()) {
  nChains <-  parallel::detectCores()
}

# set brms controls to prevent divergent transitions & improve convergence
adapt_delta <- .95
max_treedepth <- 10

###############################################################################!
# 1) Read: Raw Data ------------------------------------------------------------
###############################################################################!

#### Load data from Loaiza & Souza (2018)
# list all files for the working memory task that need to be loaded
temp <- list.files(pattern = "*RCB_WM_Aging.dat", path = here("data","RawData_LS2018"))
temp <- here("data","RawData_LS2018",temp)

# load all files and bind them by row into long-format
wmData <- do.call("rbind", lapply(temp, function(x) read.table(x, header = F)))

# name the columns contained in the data
names(wmData) <- c("id", "ageGroup", "age", "gender","test","cueCondition", "nCues", "baselineCond", "retention", "trial",
                   "item1", "item2", "item3", "item4", "item5", "colorWheelRotation",
                   "cuedItem1", "cuedItem2","targetPosition", 
                   "response", "targetColor", "NT1", "NT2", "NT3", "NT4",
                   "deviation", "recallError", "RT")

# remove practice trials
wmData <- subset(wmData, test != 0)

# remove participants with problems during data collection 
# 17 = program crash; replaced by another participant
# please see the original publication for further explanation
wmData <- subset(wmData, !(id %in% c(17)))

# convert the relevant predictors to factors and calculate the response error 
# (i.e. deviation) into radians
wmData <- wmData %>% 
  mutate(ageGroup = as.factor(ageGroup),
         RI = as.factor(retention),
         cueCond = as.factor(cueCondition),
         dev_rad = deviation/180 * pi)

# specify the levels of the experimental facctors
levels(wmData$ageGroup) <- c("Young","Old")
levels(wmData$RI) <- c("short","long")
levels(wmData$cueCond) <- c("No","Retro")

# specfiy contrasts that give equal prior width for all condition means
# for details, see: https://easystats.github.io/bayestestR/reference/contr.equalprior.html
contrasts(wmData$ageGroup) <- contr.equalprior
contrasts(wmData$RI) <- contr.equalprior
contrasts(wmData$cueCond) <- contr.equalprior

###############################################################################!
# 2) Model estimation ----------------------------------------------------------
###############################################################################!

# set up the bmmodel object
LS_model <- mixture2p(resp_error = "dev_rad")

#' set up the formula for the mixture model.
#' Although, we do not want to estimate the mean of the two von Mises distributions,
#' we have to initialize the formula to specify the dependent variable.
#' Using priors (see step 3), we will constrain the means of both von Mises distributions
#' to zero. Additionally, we will use priors to fix the precision (kappa) of the
#' second von Mises to be zero (at least practically zero). 
LS_formula <- bmf(
  # estimating fixed intercept & random intercept for kappa of the first von Mises
  kappa ~ 0 + ageGroup + ageGroup:RI + ageGroup:cueCond + ageGroup:RI:cueCond  + 
    (0 + RI + cueCond + RI:cueCond || gr(id, by = ageGroup)), 
  
  # estimating fixed intercept & random intercept for the mixing proportion 
  # of the target vonMises (i.e., p_mem)
  thetat ~ 0 + ageGroup + ageGroup:RI + ageGroup:cueCond + ageGroup:RI:cueCond  + 
    (0 + RI + cueCond + RI:cueCond || gr(id, by = ageGroup))
)

# check the default priors
default_prior(LS_formula, data = wmData, model = LS_model)

my_priors <- prior(normal(0,1), class = b, nlpar = thetat)

# fit mixture model if there is not already a results file stored
LS2018_fit <- bmm(
  # include model information
  formula = LS_formula, # specify formula for mixture model
  data    = wmData, # specify data used to estimate the mixture model
  model = LS_model,
  
  # save settings
  prior = my_priors,
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
  
  # save the results
  file = here("output","fit_E2_LS2018_bmm")
)

###############################################################################!
# 3) Model evaluation ----------------------------------------------------------
###############################################################################!


## 3.1) fit & summary ----------------------------------------------------------
# plot the posterior predictive check to evaluate overall model fit
brms::pp_check(LS2018_fit, group = "RI")

# print results summary
summary(LS2018_fit)

# test hypothesis
cue_hypothesis <- c(
  kp_cueFX_younger = "kappa_ageGroupYoung:cueCond1 = 0",
  kp_cueFX_older =  "kappa_ageGroupOld:cueCond1 = 0",
  pM_cueFX_younger = "thetat_ageGroupYoung:cueCond1 = 0",
  pM_cueFX_older = "thetat_ageGroupOld:cueCond1 = 0"
)

hypothesis(LS2018_fit, cue_hypothesis)

hyp_cue <- hypothesis(LS2018_fit, cue_hypothesis)
plot(hyp_cue)
1/hyp_cue$hypothesis$Evid.Ratio

age_hypothesis <- c(
  kp_age = "kappa_ageGroupYoung = kappa_ageGroupOld",
  pM_age = "thetat_ageGroupYoung = thetat_ageGroupOld"
)
hyp_age <- hypothesis(LS2018_fit, age_hypothesis)
1/hyp_age$hypothesis$Evid.Ratio

## 3.2) extract parameter estimates --------------------------------------------

# plot conditional effects
conditional_effects(LS2018_fit, dpar = "theta1")
conditional_effects(LS2018_fit, dpar = "kappa1")

# extract the fixed effects from the model
fixedEff <- fixef(LS2018_fit)

# determine the rows that contain the relevant parameter estimates
theta_cols <- grepl("thetat",rownames(fixedEff))
kappa_cols <- grepl("kappa",rownames(fixedEff))

# extract kappa estimates
kappa_fixedFX <- fixedEff[kappa_cols,]

# convert kappa estimates to absolute scale (radians)
kappa_fixedFX <- exp(kappa_fixedFX)

# extract theta estimates
theta_fixedFX <- fixedEff[theta_cols,]

# convert theta estimates into pMem estimates
p_Mem_fixedFX <- gtools::inv.logit(theta_fixedFX)

# print out parameter estimates
kappa_fixedFX
p_Mem_fixedFX


## 3.3) plot parameter estimates -----------------------------------------------
results_LS_2018 <- read.table(here("data","LS2018_2P_hierarchicalfit.txt"),
                              header = T, sep = ",") %>% 
  filter(param != "catActive") %>% 
  mutate(RI = retention,
         nCues = case_when(cueCond == "NoCue" ~ 0,
                           cueCond == "RetroCue" & RI == "short" ~ 1,
                           cueCond == "RetroCue" & RI == "long" ~ 2),
         ageGroup = case_when(BP_Group == "Old" ~ "old",
                              TRUE ~ "young"))



# extract posterior draws for fixed effects on kappa & theta
fixedFX_draws <- LS2018_fit %>% 
  tidy_draws() %>%
  select(starts_with("b_"),.chain,.iteration,.draw) %>% 
  mutate(kappa_old_short_no = b_kappa_ageGroupOld + 
           contrasts(wmData$RI)["short",] * `b_kappa_ageGroupOld:RI1` +
           contrasts(wmData$cueCond)["No",] * `b_kappa_ageGroupOld:cueCond1` +
           contrasts(wmData$RI)["short",] * contrasts(wmData$cueCond)["No",] * `b_kappa_ageGroupOld:RI1:cueCond1`,
         kappa_old_short_retro = b_kappa_ageGroupOld + 
           contrasts(wmData$RI)["short",] * `b_kappa_ageGroupOld:RI1` +
           contrasts(wmData$cueCond)["Retro",] * `b_kappa_ageGroupOld:cueCond1` +
           contrasts(wmData$RI)["short",] * contrasts(wmData$cueCond)["Retro",] * `b_kappa_ageGroupOld:RI1:cueCond1`,
         kappa_old_long_no = b_kappa_ageGroupOld + 
           contrasts(wmData$RI)["long",] * `b_kappa_ageGroupOld:RI1` +
           contrasts(wmData$cueCond)["No",] * `b_kappa_ageGroupOld:cueCond1` +
           contrasts(wmData$RI)["long",] * contrasts(wmData$cueCond)["No",] * `b_kappa_ageGroupOld:RI1:cueCond1`,
         kappa_old_long_retro = b_kappa_ageGroupOld + 
           contrasts(wmData$RI)["long",] * `b_kappa_ageGroupOld:RI1` +
           contrasts(wmData$cueCond)["Retro",] * `b_kappa_ageGroupOld:cueCond1` +
           contrasts(wmData$RI)["long",] * contrasts(wmData$cueCond)["Retro",] * `b_kappa_ageGroupOld:RI1:cueCond1`,
         kappa_young_short_no = b_kappa_ageGroupYoung + 
           contrasts(wmData$RI)["short",] * `b_kappa_ageGroupYoung:RI1` +
           contrasts(wmData$cueCond)["No",] * `b_kappa_ageGroupYoung:cueCond1` +
           contrasts(wmData$RI)["short",] * contrasts(wmData$cueCond)["No",] * `b_kappa_ageGroupYoung:RI1:cueCond1`,
         kappa_young_short_retro = b_kappa_ageGroupYoung + 
           contrasts(wmData$RI)["short",] * `b_kappa_ageGroupYoung:RI1` +
           contrasts(wmData$cueCond)["Retro",] * `b_kappa_ageGroupYoung:cueCond1` +
           contrasts(wmData$RI)["short",] * contrasts(wmData$cueCond)["Retro",] * `b_kappa_ageGroupYoung:RI1:cueCond1`,
         kappa_young_long_no = b_kappa_ageGroupYoung + 
           contrasts(wmData$RI)["long",] * `b_kappa_ageGroupYoung:RI1` +
           contrasts(wmData$cueCond)["No",] * `b_kappa_ageGroupYoung:cueCond1` +
           contrasts(wmData$RI)["long",] * contrasts(wmData$cueCond)["No",] * `b_kappa_ageGroupYoung:RI1:cueCond1`,
         kappa_young_long_retro = b_kappa_ageGroupYoung + 
           contrasts(wmData$RI)["long",] * `b_kappa_ageGroupYoung:RI1` +
           contrasts(wmData$cueCond)["Retro",] * `b_kappa_ageGroupYoung:cueCond1` +
           contrasts(wmData$RI)["long",] * contrasts(wmData$cueCond)["Retro",] * `b_kappa_ageGroupYoung:RI1:cueCond1`,
         thetat_old_short_no = b_thetat_ageGroupOld + 
           contrasts(wmData$RI)["short",] * `b_thetat_ageGroupOld:RI1` +
           contrasts(wmData$cueCond)["No",] * `b_thetat_ageGroupOld:cueCond1` +
           contrasts(wmData$RI)["short",] * contrasts(wmData$cueCond)["No",] * `b_thetat_ageGroupOld:RI1:cueCond1`,
         thetat_old_short_retro = b_thetat_ageGroupOld + 
           contrasts(wmData$RI)["short",] * `b_thetat_ageGroupOld:RI1` +
           contrasts(wmData$cueCond)["Retro",] * `b_thetat_ageGroupOld:cueCond1` +
           contrasts(wmData$RI)["short",] * contrasts(wmData$cueCond)["Retro",] * `b_thetat_ageGroupOld:RI1:cueCond1`,
         thetat_old_long_no = b_thetat_ageGroupOld + 
           contrasts(wmData$RI)["long",] * `b_thetat_ageGroupOld:RI1` +
           contrasts(wmData$cueCond)["No",] * `b_thetat_ageGroupOld:cueCond1` +
           contrasts(wmData$RI)["long",] * contrasts(wmData$cueCond)["No",] * `b_thetat_ageGroupOld:RI1:cueCond1`,
         thetat_old_long_retro = b_thetat_ageGroupOld + 
           contrasts(wmData$RI)["long",] * `b_thetat_ageGroupOld:RI1` +
           contrasts(wmData$cueCond)["Retro",] * `b_thetat_ageGroupOld:cueCond1` +
           contrasts(wmData$RI)["long",] * contrasts(wmData$cueCond)["Retro",] * `b_thetat_ageGroupOld:RI1:cueCond1`,
         thetat_young_short_no = b_thetat_ageGroupYoung + 
           contrasts(wmData$RI)["short",] * `b_thetat_ageGroupYoung:RI1` +
           contrasts(wmData$cueCond)["No",] * `b_thetat_ageGroupYoung:cueCond1` +
           contrasts(wmData$RI)["short",] * contrasts(wmData$cueCond)["No",] * `b_thetat_ageGroupYoung:RI1:cueCond1`,
         thetat_young_short_retro = b_thetat_ageGroupYoung + 
           contrasts(wmData$RI)["short",] * `b_thetat_ageGroupYoung:RI1` +
           contrasts(wmData$cueCond)["Retro",] * `b_thetat_ageGroupYoung:cueCond1` +
           contrasts(wmData$RI)["short",] * contrasts(wmData$cueCond)["Retro",] * `b_thetat_ageGroupYoung:RI1:cueCond1`,
         thetat_young_long_no = b_thetat_ageGroupYoung + 
           contrasts(wmData$RI)["long",] * `b_thetat_ageGroupYoung:RI1` +
           contrasts(wmData$cueCond)["No",] * `b_thetat_ageGroupYoung:cueCond1` +
           contrasts(wmData$RI)["long",] * contrasts(wmData$cueCond)["No",] * `b_thetat_ageGroupYoung:RI1:cueCond1`,
         thetat_young_long_retro = b_thetat_ageGroupYoung + 
           contrasts(wmData$RI)["long",] * `b_thetat_ageGroupYoung:RI1` +
           contrasts(wmData$cueCond)["Retro",] * `b_thetat_ageGroupYoung:cueCond1` +
           contrasts(wmData$RI)["long",] * contrasts(wmData$cueCond)["Retro",] * `b_thetat_ageGroupYoung:RI1:cueCond1`) %>% 
  select(.chain,.iteration,.draw, starts_with("kappa"), starts_with("thetat")) %>% 
  pivot_longer(cols = c(starts_with("kappa"), starts_with("thetat")),
               names_to = "modelPar",
               values_to = "postSample") %>% 
  mutate(par = str_split_i(modelPar,"_",1),
         ageGroup = str_split_i(modelPar,"_",2),
         RI = str_split_i(modelPar,"_",3),
         RetroCue = str_split_i(modelPar,"_",4)) %>% 
  select(-modelPar) %>% 
  filter(par == "kappa" | par == "thetat") %>% 
  mutate(postSample_abs = case_when(par == "kappa" ~ (sqrt(1/exp(postSample))/pi) * 180,
                                    par == "thetat" ~ inv_logit_scaled(postSample)),
         nCues = case_when(RetroCue == "no" ~ 0,
                           RetroCue == "retro" & RI == "short" ~ 1,
                           RetroCue == "retro" & RI == "long" ~ 2))

# plot kappa results
kappa_plot <- ggplot(data = fixedFX_draws %>% filter(par == "kappa"),
                     aes(x = RI, y = postSample_abs, 
                         color = as.factor(nCues), shape = as.factor(nCues))) +
  facet_grid(. ~ ageGroup) +
  coord_cartesian(ylim = c(10,35)) +
  geom_half_violin(position = position_nudge(x = .1, y = 0), aes(fill = as.factor(nCues)), side = "r",
                   adjust = 1, trim = TRUE, alpha = 0.9, show.legend = FALSE, scale = "width") +
  stat_summary(geom = "pointrange", fun.data = mean_hdci,
               size = 0.7, linewidth = 0.8,
               position = position_dodge(0.1)) +
  geom_point(data = results_LS_2018 %>% filter(param == "contSD"), 
             aes(y = mean, x = RI, color = as.factor(nCues)),
             shape = "diamond", size = 3,
             position = position_nudge(x = -.1, y = 0)) +
  scale_fill_grey(start = 0, end = .7) +
  scale_color_grey(start = 0, end = .7) +
  labs(x = "Retention Interval", y = "Memory imprecision (SD)", 
       fill = "No. of Cues", color = "No. of Cues", shape = "No. of Cues",
       title = "B") +
  guides(color = "none", shape = "none") +
  clean_plot()
kappa_plot

# plot pMem results
pMem_plot <- ggplot(data = fixedFX_draws %>% filter(par == "thetat"),
                    aes(x = RI, y = postSample_abs, 
                        color = as.factor(nCues), shape = as.factor(nCues))) +
  facet_grid(. ~ ageGroup) +
  theme(legend.position = c(0.2, 0.8)) +
  coord_cartesian(ylim = c(0.35,1)) +
  geom_half_violin(position = position_nudge(x = .1, y = 0), aes(fill = as.factor(nCues)), side = "r",
                   adjust = 1.2, trim = FALSE, alpha = 0.9, colour = NA, show.legend = FALSE, scale = "width") +
  stat_summary(geom = "pointrange", fun.data = mean_hdci,
               size = 0.7, linewidth = 0.8,
               position = position_dodge(0.1)) +
  geom_point(data = results_LS_2018 %>% filter(param == "pMem"), 
             aes(y = mean, x = RI, color = as.factor(nCues)),
             shape = "diamond", size = 3,
             position = position_nudge(x = -.1, y = 0)) +
  scale_fill_grey(start = 0, end = 0.7) +
  scale_color_grey(start = 0, end = 0.7) +
  labs(x = "Retention Interval", y = expression(P[mem]), 
       color = "No. of Cues", shape = "No. of Cues",
       title = "A") +
  clean_plot()
pMem_plot

# patch plots together
joint_plot <- (pMem_plot | kappa_plot)

# show joint plot
joint_plot

# save plots with high resolution
ggsave(
  filename = here("figures/plot_kappaEst_LS2018.jpeg"),
  plot = kappa_plot, width = 6, height = 6
)
ggsave(
  filename = here("figures/plot_pmemEst_LS2018.jpeg"),
  plot = pMem_plot, width = 6, height = 6
)

ggsave(
  filename = here("figures/plot_jointRes_LS2018.jpeg"),
  plot = joint_plot, width = 5*2, height = 5
)
