# start fresh
rm(list = ls()) # clean up work space
graphics.off()  # switch off graphics device

# load required packages
pacman::p_load(here, osfr)

# retrieve OSF Project of the Tutorial Paper
Tutorial_project <- osf_retrieve_node("vsrz4")

# Download all missing files to the output folder
Tutorial_project %>%
  osf_ls_files(n_max = 20) %>% 
  osf_download(path = here("output"), conflicts = "skip", verbose = FALSE)

# Print short message about the files being downloaded
print("All missing results files have been downloaded sucesfully!")