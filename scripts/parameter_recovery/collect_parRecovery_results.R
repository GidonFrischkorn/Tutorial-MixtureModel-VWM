library(here)
library(glue)
library(tidytable)

# list all recovery files
recovery_files <- list.files(here("output","recovery_results"), pattern = "par_rec_fits", full.names = TRUE)

# read all recovery files
fits <- lapply(recovery_files, readRDS)

nReplications <- length(fits)

df_hyperPar_gen <- df_hyperPar_ml <- df_hyperPar_bmm <- data.table(
  pmem = numeric(),
  kappa = numeric(),
  n_trials = numeric(),
  n_sub = numeric(),
  n_rep = integer()
)

df_subPar_gen <- df_subPar_ml <- df_subPar_bmm <- data.table(
  id = integer(),
  pmem = numeric(),
  kappa = numeric(),
  n_trials = numeric(),
  n_sub = numeric(),
  n_rep = integer()
)


for (i in 1:nReplications) {
  cond_hyperPar_gen <- fits[[i]]$hyperPars %>% t() %>% data.table(.) %>% mutate(n_rep = i)
  cond_hyperPar_gen$kappa <- exp(cond_hyperPar_gen$kappa)
  cond_hyperPar_gen$pmem <- exp(cond_hyperPar_gen$pmem)/(1+exp(cond_hyperPar_gen$pmem))
  df_hyperPar_gen <- rbind(df_hyperPar_gen, cond_hyperPar_gen)
  
  cond_subPar_gen <- fits[[i]]$pars %>% mutate(n_rep = i)
  df_subPar_gen <- rbind(df_subPar_gen, cond_subPar_gen)
  
  cond_subPar_ml <- fits[[i]]$par_ml %>% mutate(n_rep = i)
  cond_subPar_ml$id <- as.numeric(cond_subPar_ml$id)
  df_subPar_ml <- rbind(df_subPar_ml, cond_subPar_ml)
  
  cond_subPar_bmm <- fits[[i]]$par_bmm %>% mutate(n_rep = i)
  df_subPar_bmm <- rbind(df_subPar_bmm, cond_subPar_bmm)
  
  cond_hyperPar_bmm <- fits[[i]]$hyperPar_bmm %>% data.table(.) %>% mutate(n_rep = i)
  colnames(cond_hyperPar_bmm) <- c("kappa","pmem", "n_sub", "n_trials","n_rep")
  cond_hyperPar_bmm$kappa <- exp(cond_hyperPar_bmm$kappa)
  cond_hyperPar_bmm$pmem <- exp(cond_hyperPar_bmm$pmem)/(1+exp(cond_hyperPar_bmm$pmem))
  df_hyperPar_bmm <- rbind(df_hyperPar_bmm, cond_hyperPar_bmm)
}
rm(cond_hyperPar_bmm, cond_hyperPar_gen, cond_subPar_gen, cond_subPar_ml, cond_subPar_bmm)

df_hyperPar_ml <- df_subPar_ml %>% 
  summarise(kappa = mean(kappa),
            pmem = mean(pmem),
            .by = c("n_rep", "n_trials","n_sub"))

df_hyperPar_bmm <- df_hyperPar_bmm %>% 
  pivot_longer(cols = c("pmem", "kappa"), names_to = "parameter", values_to = "est_bmm")
df_hyperPar_ml <- df_hyperPar_ml %>%
  pivot_longer(cols = c("pmem", "kappa"), names_to = "parameter", values_to = "est_ml")
df_hyperPar_gen <- df_hyperPar_gen %>%
  pivot_longer(cols = c("pmem", "kappa"), names_to = "parameter", values_to = "gen")
df_hyperPar <- df_hyperPar_gen %>% 
  left_join(df_hyperPar_ml, by = c("n_rep", "n_sub", "n_trials", "parameter")) %>% 
  left_join(df_hyperPar_bmm, by = c("n_rep", "n_sub", "n_trials", "parameter"))
rm(df_hyperPar_bmm, df_hyperPar_ml, df_hyperPar_gen)

df_subPar_gen <- df_subPar_gen %>% 
  pivot_longer(cols = c("pmem", "kappa"), names_to = "parameter", values_to = "gen")
df_subPar_ml <- df_subPar_ml %>%
  pivot_longer(cols = c("pmem", "kappa"), names_to = "parameter", values_to = "est_ml")
df_subPar_bmm <- df_subPar_bmm %>%
  pivot_longer(cols = c("pmem", "kappa"), names_to = "parameter", values_to = "est_bmm")
df_subPar <- df_subPar_gen %>%
  left_join(df_subPar_ml, by = c("id", "n_rep", "n_sub", "n_trials", "parameter")) %>% 
  left_join(df_subPar_bmm, by = c("id", "n_rep", "n_sub", "n_trials", "parameter"))
rm(df_subPar_gen, df_subPar_ml, df_subPar_bmm)
