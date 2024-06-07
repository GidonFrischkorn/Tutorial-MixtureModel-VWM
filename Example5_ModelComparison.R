#' This is the tutorial script for setting up the Interference Measurement Model (IMMfull)
#' for visual working memory tasks that use continuous report recall procedures.
#' 
#' In this script, you will see:
#'  1) how the model is set up using the bmm package, 
#'  2) how a simple version of the model is estimates, and 
#'  3) how the model can be evaluated and results extracted and plotted.

# start fresh
rm(list = ls()) # clean up work space
graphics.off()  # switch off graphics device

# load required packages
pacman::p_load(here, bmm, brms, tidyverse, tidybayes, patchwork, gghalves)


# 3pMM_fit <- readRDS(here("output","fit_E5_OL2017_3pMM.rds"))
# imm_bsc_fit <- readRDS(here("output","fit_E5_OL2017_IMMbsc.rds"))
# imm_abc_fit <- readRDS(here("output","fit_E5_OL2017_IMMabc.rds")) 
# imm_full_fit <- readRDS(here("output","fit_E5_OL2017_IMMfull.rds"))
# 
# bridge_imm_full <- bridge_sampler(imm_full_fit, repetition = 20, cores = 6)
# bridge_imm_bsc <- bridge_sampler(imm_bsc_fit, repetition = 20, cores = 6)
# bayes_factor(bridge_imm_full, bridge_imm_bsc)

mixture3p_fit <- readRDS(file = here("output","fit_E5_OL2017_3pMM.rds"))

if (!file.exists(here("output","E5_bridge_3p.rds"))) {
  bridge_3p <- bridge_sampler(mixture3p_fit, repetitions = 10, cores = 4)
  saveRDS(bridge_3p, file = here("output","E5_bridge_3p.rds"))
} else {
  bridge_3p <- readRDS(here("output","E5_bridge_3p.rds"))
}

if (!file.exists(here("output","E5_bridge_IMMabc.rds"))) {
  bridge_imm_abc <- bridge_sampler(imm_abc_fit, repetitions = 10, cores = 4)
  saveRDS(bridge_imm_abc, file = here("output","E5_bridge_IMMabc.rds"))
} else {
  bridge_imm_abc <- readRDS(here("output","E5_bridge_IMMabc.rds"))
}

bf <- bayes_factor(bridge_3p, bridge_imm_abc)
hist(bf$bf)