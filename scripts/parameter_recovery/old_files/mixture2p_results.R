###############################################################################-
#### SETUP ####
###############################################################################-

extract_fit_stats <- function(fits) {
  stats <- lapply(fits, par_rec_stats)
  stats <- do.call(rbind, stats)
  stats$replicate <- rep(1:100, each = 4)
  stats
}

library(bmm)
library(tidyverse)
library(ggdist)
library(here)
library(patchwork)

source(here("scripts/parameter_recovery/mixture2p_functions.R"))

fits50 <- readRDS(here("output/par_rec_fits50.rds"))
fits100 <- readRDS(here("output/par_rec_fits100.rds"))





###############################################################################-
#### RESULTS ####
###############################################################################-

stats50 <- extract_fit_stats(fits50)

(plot_r_50 <- stats50 %>% 
  ggplot(aes(r, param, fill = engine)) +
  stat_halfeye(position = "dodge", scale = 0.5) +
  stat_dotsinterval(position = "dodge", scale = 0.5, side = "bottom") +
  theme_ggdist() +
  facet_wrap(~param, scale = 'free') +
  xlab('Correlation between subject-level\ntrue and estimated parameters') +
  scale_y_discrete("", labels = "") +
  ggtitle('50 observations per participant\n\nA) Correlation'))


(plot_rmse_50 <- stats50 %>% 
  ggplot(aes(rmse, param, fill = engine)) +
  stat_halfeye(position = "dodge", scale = 0.5, height = 4, normalize = "groups") +
  stat_dotsinterval(position = "dodge", scale = 0.5, height = 4, 
                    side = "bottom", normalize = "groups") +
  theme_ggdist() +
  facet_wrap(~param, scale = 'free') +
  xlab('Root mean squared error between\nsubject-level true and estimated parameters') +
  scale_y_discrete("", labels = "") +
  ggtitle('\n\nB) RMSE'))

(plot_pop_mean_50 <- stats50 %>% 
  ggplot(aes(pop_par_diff, param, fill = engine)) +
  stat_halfeye(position = "dodge", scale = 0.5, height = 4, normalize = "groups") +
  stat_dotsinterval(position = "dodge", scale = 0.5, height = 4, 
                    side = "bottom", normalize = "groups") +
  theme_ggdist() +
  facet_wrap(~param, scale = 'free') +
  xlab('Error in estimating\nthe population mean parameter') +
  scale_y_discrete("", labels = "") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  ggtitle("\n\nC) Population mean parameter estimation error"))

plot_r_50 + theme(legend.position = "none") + 
  plot_rmse_50 + theme(legend.position = "none") +
  plot_pop_mean_50


stats100 <- extract_fit_stats(fits100)

(plot_r_100 <- stats100 %>% 
    ggplot(aes(r, param, fill = engine)) +
    stat_halfeye(position = "dodge", scale = 0.5) +
    stat_dotsinterval(position = "dodge", scale = 0.5, side = "bottom") +
    theme_ggdist() +
    facet_wrap(~param, scale = 'free') +
    xlab('Correlation between subject-level\ntrue and estimated parameters') +
    scale_y_discrete("", labels = "") +
    ggtitle('\n100 observations per participant\n\nD) Correlation'))


(plot_rmse_100 <- stats100 %>% 
    ggplot(aes(rmse, param, fill = engine)) +
    stat_halfeye(position = "dodge", scale = 0.5, height = 4, normalize = "groups") +
    stat_dotsinterval(position = "dodge", scale = 0.5, height = 4, 
                      side = "bottom", normalize = "groups") +
    theme_ggdist() +
    facet_wrap(~param, scale = 'free') +
    xlab('Root mean squared error between\nsubject-level true and estimated parameters') +
    scale_y_discrete("", labels = "") +
    ggtitle('\n\n\nE) RMSE'))

(plot_pop_mean_100 <- stats100 %>% 
    ggplot(aes(pop_par_diff, param, fill = engine)) +
    stat_halfeye(position = "dodge", scale = 0.5, height = 4, normalize = "groups") +
    stat_dotsinterval(position = "dodge", scale = 0.5, height = 4, 
                      side = "bottom", normalize = "groups") +
    theme_ggdist() +
    facet_wrap(~param, scale = 'free') +
    xlab('Error in estimating\nthe population mean parameter') +
    scale_y_discrete("", labels = "") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    ggtitle("\n\n\nF) Population mean parameter estimation error"))

plot_r_100 + theme(legend.position = "none") + 
  plot_rmse_100 + theme(legend.position = "none") +
  plot_pop_mean_100



(plot_r_50 + theme(legend.position = "none") + 
  plot_rmse_50 + theme(legend.position = "none") +
  plot_pop_mean_50) /
  (plot_r_100 + theme(legend.position = "none") + 
   plot_rmse_100 + theme(legend.position = "none") +
   plot_pop_mean_100)


###############################################################################-
#### DIFFERENCES IN CORRELATION ####
###############################################################################-

corr_diff_50 <- stats50 %>% 
  select(engine, param, r) %>% 
  pivot_wider(names_from = engine, values_from = r) %>% 
  unnest() %>% 
  mutate(diff = BMM - ML)

(plot_corr_diff_50 <- corr_diff_50 %>%
  ggplot(aes(diff, param)) +
  stat_halfeye(position = "dodge", scale = 0.5) +
  stat_dotsinterval(position = "dodge", scale = 0.5, side = "bottom") +
  theme_ggdist() +
  facet_wrap(~param, scale = 'free') +
  xlab('Difference in correlation between\nsubject-level true and estimated parameters') +
  scale_y_discrete("", labels = "") +
  ggtitle('50 observations per participant\n\n Correlation difference'))

corr_diff_50 %>% 
  group_by(param) %>% 
  reframe(rdiff = quantile(diff, c(0.025, 0.250, 0.5, 0.75, 0.975))) %>%
  ungroup() %>% 
  mutate(engine = "BMM - ML",
         quantile = rep(c("Q2.5", "Q25", "Q50", "Q75", "Q97.5"), nrow(.)/5)) %>%
  pivot_wider(names_from = quantile, values_from = rdiff)
