# Reproducing results from Loaiza & Souza (2018)
# "Is Refreshing in Working Memory Impaired in Older Age? Evidence from the Retro-Cue Paradigm"

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

###############################################################################!
# 2) Model estimation ----------------------------------------------------------
###############################################################################!

#' First we specify a mixture of von Mises distributions.
#' The first distribution is the memory distribution and
#' the second distribution is the guessing
LS_mixFamily <- mixture(von_mises,von_mises, order =  "none")

#' Then, we set up the formula for the mixture model.
#' Although, we do not want to estimate the mean of the two von Mises distributions,
#' we have to initialize the formula to specify the dependent variable.
#' Using priors (see step 3), we will constrain the means of both von Mises distributions
#' to zero. Additionally, we will use priors to fix the precision (kappa) of the
#' second von Mises to be zero (at least practically zero). 
LS_mixFormula <- bf(dev_rad ~ 1,    # Initializing the dependent variable
                    mu2 ~ 1,
                    # estimating fixed intercept & random intercept for kappa of the first von Mises
                    kappa1 ~ 0 + ageGroup:RI:cueCond + (0 + RI:cueCond || gr(id, by = ageGroup)), 
                    kappa2 ~ 1,
                    # estimating fixed intercept & random intercept for the mixing proportion 
                    # of the first vonMises (i.e., p_mem)
                    theta1 ~ 0 + ageGroup:RI:cueCond + (0 + RI:cueCond || gr(id, by = ageGroup)))

# constrain parameters using priors
LS_mixPriors <- 
  # fix mean of the first von Mises to zero
  prior(constant(0), class = Intercept, dpar = "mu1") +
  # fix mean of the second von Mises to zero
  prior(constant(0), class = Intercept, dpar = "mu2") +
  # fix kappa of the second von Mises to (alomst) zero
  prior(constant(-100), class = Intercept, dpar = "kappa2") +
  prior(normal(0,0.5), class = "b", dpar = "theta1") +
  prior(normal(0,0.5), class = "b", dpar = "kappa1")

# fit mixture model if there is not already a results file stored
if (!file.exists(here("output","fit_E2_LS2018.RData"))) {
  # fit the mixture model using brms
  fit_LS2018_mixModel <- brm(
    # include model information
    formula = LS_mixFormula, # specify formula for mixture model
    data    = wmData, # specify data used to estimate the mixture model
    family  = LS_mixFamily, # call the defined mixture family
    prior   = LS_mixPriors, # use the used defined priors,
  
    # save settings
    sample_prior = TRUE,
    save_pars = save_pars(all = TRUE),
    
    # add brms settings
    warmup = warmup_samples,
    iter = warmup_samples + postwarmup_samples, 
    chains = nChains,
    
    # control commands for the sampler
    control = list(adapt_delta = adapt_delta, 
                   max_treedepth = max_treedepth)
  )
  
  # save results into file
  save(fit_LS2018_mixModel, 
       file = here("output","fit_E2_LS2018.RData"),
       compress = "xz")
  
} else {
  # load results file
  load(file = here("output","fit_E2_LS2018.RData"))
}

###############################################################################!
# 3) Model evaluation ----------------------------------------------------------
###############################################################################!

## 3.1) fit & summary ----------------------------------------------------------
# plot the posterior predictive check to evaluate overall model fit
pp_check(fit_LS2018_mixModel)

# print results summary
summary(fit_LS2018_mixModel)

# test hypothesis
hypothesis(fit_LS2018_mixModel,
           c(hyp1 = "theta1_ageGroupYoung:RIshort:cueCondNo > theta1_ageGroupOld:RIshort:cueCondNo"))

## 3.2) extract parameter estimates --------------------------------------------

# extract the fixed effects from the model
fixedEff <- fixef(fit_LS2018_mixModel)

# determine the rows that contain the relevant parameter estimates
theta_cols <- grepl("theta",rownames(fixedEff))
kappa_cols <- grepl("kappa1",rownames(fixedEff))

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
         ageGroup = case_when(BP_Group == "Old" ~ "Old",
                              TRUE ~ "Young"))
  
# extract posterior draws for fixed effects on kappa & theta
fixedFX_draws <- fit_LS2018_mixModel %>% 
  tidy_draws() %>%
  select(starts_with("b_"),.chain,.iteration,.draw) %>% 
  pivot_longer(cols = starts_with("b_"),
               names_to = "modelPar",
               values_to = "postSample") %>% 
  mutate(par = str_split_i(modelPar,"_",2),
         conds = str_split_i(modelPar,"_",3),
         ageGroup = str_split_i(conds,":",1),
         RI = str_split_i(conds,":",2),
         RetroCue = str_split_i(conds,":",3),
         ageGroup = str_remove(ageGroup,"ageGroup"),
         RI = str_remove(RI,"RI"),
         RetroCue = str_remove(RetroCue,"cueCond")) %>% 
  select(-modelPar, -conds) %>% 
  filter(par == "kappa1" | par == "theta1") %>% 
  mutate(postSample_abs = case_when(par == "kappa1" ~ (sqrt(1/exp(postSample))/pi) * 180,
                                    par == "theta1" ~ inv_logit_scaled(postSample)),
         nCues = case_when(RetroCue == "No" ~ 0,
                           RetroCue == "Retro" & RI == "short" ~ 1,
                           RetroCue == "Retro" & RI == "long" ~ 2))

# plot kappa results
kappa_plot <- ggplot(data = fixedFX_draws %>% filter(par == "kappa1"),
                     aes(x = RI, y = postSample_abs, color = as.factor(nCues))) +
  facet_grid(. ~ ageGroup) +
  coord_cartesian(ylim = c(5,50)) +
  geom_half_violin(position = position_nudge(x = .1, y = 0), aes(fill = as.factor(nCues)), side = "r",
                   adjust = 1, trim = TRUE, alpha = 0.9, show.legend = FALSE, scale = "width") +
  stat_summary(geom = "pointrange", fun.data = mode_hdi,
               size = 0.3, linewidth = 0.8,
               position = position_dodge(0.1)) +
  geom_point(data = results_LS_2018 %>% filter(param == "contSD"), 
             aes(y = mean, x = RI, color = as.factor(nCues)),
             shape = "diamond", size = 2.5,
             position = position_nudge(x = -.1, y = 0)) +
  scale_fill_grey(start = 0, end = .8) +
  scale_color_grey(start = 0, end = .8) +
  labs(x = "Retention Interval", y = "Memory imprecision (SD)", fill = "No. of Cues", color = "No. of Cues",
       title = "B") +
  guides(color = "none") +
  clean_plot()
kappa_plot

# plot pMem results
pMem_plot <- ggplot(data = fixedFX_draws %>% filter(par == "theta1"),
                     aes(x = RI, y = postSample_abs, color = as.factor(nCues))) +
  facet_grid(. ~ ageGroup) +
  theme(legend.position = c(0.25, 0.8)) +
  coord_cartesian(ylim = c(0.35,1)) +
  geom_half_violin(position = position_nudge(x = .1, y = 0), aes(fill = as.factor(nCues)), side = "r",
                   adjust = 1.2, trim = FALSE, alpha = 0.9, colour = NA, show.legend = FALSE, scale = "width") +
  stat_summary(geom = "pointrange", fun.data = mean_hdi,
               size = 0.3, linewidth = 0.8,
               position = position_dodge(0.1)) +
  geom_point(data = results_LS_2018 %>% filter(param == "pMem"), 
             aes(y = mean, x = RI, color = as.factor(nCues)),
             shape = "diamond", size = 2.5,
             position = position_nudge(x = -.1, y = 0)) +
  scale_fill_grey(start = 0, end = 0.8) +
  scale_color_grey(start = 0, end = 0.8) +
  labs(x = "Retention Interval", y = expression(P[mem]), fill = "No. of Cues", color = "No. of Cues",
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
  plot = joint_plot, width = 6*2, height = 6
)
