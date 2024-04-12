# clean up work space at start up
rm(list = ls())

# switch off all graphics devices
graphics.off()

# load required packages
pacman::p_load(rtdists, brms, here, tidyverse, SimDesign)

# source all files in the functions folder
files.sources = list.files(here("functions"))
files.sources <- paste("functions",files.sources, sep = "/")
sapply(files.sources, source)

recSim_dm <- createDesign(
  nSub = c(25,50,100,200),
  nTrials = c(50,100,200),
  method = c("ml", "bmm")
)

SimFunctions(filename = "functions_recDMsim.R")

