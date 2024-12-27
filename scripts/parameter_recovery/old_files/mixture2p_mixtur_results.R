###############################################################################-
#### setup ####
###############################################################################-

library(bmm)
library(tidyverse)
library(ggdist)
library(patchwork)
library(scales)

# generating parameters
pars <- expand.grid(kappa = seq(1, 16, 0.5),
                    pmem = seq(0.6, 1, 0.025))
pars$id <- as.character(1:nrow(pars))

# TODO: read the dataset once fitting is finished (for now from tmp.RData)
load('output/mixture2p_mixtur_grid_recovery.RData')

# convert data into a better format
fits20 <- lapply(1:ncol(fits20), function(x) as.data.frame(fits20[,x]))
fits20 <- do.call(rbind, fits20)
fits20$repl <- rep(1:200, each = nrow(pars))
fits50 <- lapply(1:ncol(fits50), function(x) as.data.frame(fits50[,x]))
fits50 <- do.call(rbind, fits50)
fits50$repl <- rep(1:200, each = nrow(pars))
fits100 <- lapply(1:ncol(fits100), function(x) as.data.frame(fits100[,x]))
fits100 <- do.call(rbind, fits100)
fits100$repl <- rep(1:200, each = nrow(pars))
fits200 <- lapply(1:ncol(fits200), function(x) as.data.frame(fits200[,x]))
fits200 <- do.call(rbind, fits200)
fits200$repl <- rep(1:200, each = nrow(pars))
fits500 <- lapply(1:ncol(fits500), function(x) as.data.frame(fits500[,x]))
fits500 <- do.call(rbind, fits500)
fits500$repl <- rep(1:200, each = nrow(pars))

# add generating parms
fits20 <- left_join(fits20, pars, by = "id")
fits50 <- left_join(fits50, pars, by = "id")
fits100 <- left_join(fits100, pars, by = "id")
fits200 <- left_join(fits200, pars, by = "id")
fits500 <- left_join(fits500, pars, by = "id")

###############################################################################-
#### PLOTS ####
###############################################################################-

# plot some fits
fits20 %>% 
  ggplot(aes(kappa.y, kappa.x)) +
  geom_point() +
  facet_wrap(~pmem.y, scales = "free")

fits20 %>% 
  ggplot(aes(kappa.y, kappa.x)) +
  coord_cartesian(ylim = c(0, 50), xlim = c(0, 50)) +
  geom_point(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  facet_wrap(~pmem.y, scales = "free")

fits20 %>% 
  ggplot(aes(pmem.y, pmem.x)) +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) +
  geom_point(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  facet_wrap(~kappa.y, scales = "free")


# plot some fits
fits50 %>% 
  ggplot(aes(kappa.y, kappa.x)) +
  geom_point() +
  facet_wrap(~pmem.y, scales = "free")

fits50 %>% 
  ggplot(aes(kappa.y, kappa.x)) +
  coord_cartesian(ylim = c(0, 50), xlim = c(0, 50)) +
  geom_point(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  facet_wrap(~pmem.y, scales = "free")

fits50 %>% 
  ggplot(aes(pmem.y, pmem.x)) +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) +
  geom_point(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  facet_wrap(~kappa.y, scales = "free")


# plot some fits
fits100 %>% 
  ggplot(aes(kappa.y, kappa.x)) +
  geom_point() +
  facet_wrap(~pmem.y, scales = "free")

fits100 %>% 
  ggplot(aes(kappa.y, kappa.x)) +
  coord_cartesian(ylim = c(0, 50), xlim = c(0, 50)) +
  geom_point(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  facet_wrap(~pmem.y, scales = "free")

fits100 %>% 
  ggplot(aes(pmem.y, pmem.x)) +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) +
  geom_point(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  facet_wrap(~kappa.y, scales = "free")


# plot some fits
fits200 %>% 
  ggplot(aes(kappa.y, kappa.x)) +
  geom_point() +
  facet_wrap(~pmem.y, scales = "free")

fits200 %>% 
  ggplot(aes(kappa.y, kappa.x)) +
  coord_cartesian(ylim = c(0, 50), xlim = c(0, 50)) +
  geom_point(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  facet_wrap(~pmem.y, scales = "free")

fits200 %>% 
  ggplot(aes(pmem.y, pmem.x)) +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) +
  geom_point(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  facet_wrap(~kappa.y, scales = "free")

# plot some fits
fits500 %>% 
  ggplot(aes(kappa.y, kappa.x)) +
  geom_point() +
  facet_wrap(~pmem.y, scales = "free")

fits500 %>% 
  ggplot(aes(kappa.y, kappa.x)) +
  coord_cartesian(ylim = c(0, 50), xlim = c(0, 50)) +
  geom_point(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  facet_wrap(~pmem.y, scales = "free")

fits500 %>% 
  ggplot(aes(pmem.y, pmem.x)) +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) +
  geom_point(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  facet_wrap(~kappa.y, scales = "free")


###############################################################################-
#### STATS ####
###############################################################################-

par_rec_stats <- function(x) {
  x %>% 
    group_by(repl) %>% 
    summarise(
      r_kappa = cor(kappa.y, kappa.x),
      r_pmem = cor(pmem.y, pmem.x),
      rmse_kappa = sqrt(mean((kappa.y - kappa.x)^2)),
      rmse_pmem = sqrt(mean((pmem.y - pmem.x)^2))
    ) %>% 
    pivot_longer(cols = r_kappa:rmse_pmem, names_to = c('stat','par'), 
                 names_sep = "_", values_to = "value")
}

stats20 <- par_rec_stats(fits20)
stats50 <- par_rec_stats(fits50)
stats100 <- par_rec_stats(fits100)
stats200 <- par_rec_stats(fits200)
stats500 <- par_rec_stats(fits500)

stats20$Nobs <- 20
stats50$Nobs <- 50
stats100$Nobs <- 100
stats200$Nobs <- 200
stats500$Nobs <- 500

stats <- rbind(stats20, stats50, stats100, stats200, stats500)
stats$Nobs <- as.factor(stats$Nobs)
stats <- stats %>% 
  mutate(stat = recode(stat, r = "correlation", rmse = "RMSE"))

###############################################################################-
#### STATS ANALYSIS ####
###############################################################################-

stats %>% 
  group_by(Nobs, stat, par) %>% 
  summarise(
    mean = mean(value),
    sd = sd(value),
    median = median(value),
    HDI95 = capture.output(bayestestR::hdi(value)),
    HDI50 = capture.output(bayestestR::hdi(value, ci = 0.5)),
  ) %>% 
  arrange(stat, par, Nobs)


(stats_res <- stats %>% 
  group_by(Nobs, stat, par) %>% 
  summarise(value = posterior::rvar(value)) %>% 
  arrange(stat, par, Nobs) %>% 
  mutate(HDI95 = {
    hdi <- bayestestR::hdi(value)
    paste0("[",sprintf("%.2f", hdi$CI_low), ", ", sprintf("%.2f", hdi$CI_high), "]")
  },
  median = median(value)
  ))

###############################################################################-
#### STATS PLOTS ####
###############################################################################-

(p1 <- stats %>%
  filter(par == "kappa", stat == "correlation") %>%
  ggplot(aes(value, y = Nobs)) +
  geom_segment(
    data = filter(stats_res, par == "kappa", stat == "correlation"),
    aes(x = median, yend = -Inf),
    color = "gray65",
    linetype = "dashed"
  ) +
  stat_halfeye(normalize = "groups", fill = hue_pal()(2)[1]) +
  facet_grid(cols = vars(par),
             scales = "free_x") +
  scale_x_continuous(
    "Correlation (regular spacing)",
    breaks = filter(stats_res, par == "kappa", stat == "correlation")$median[1:4],
    labels = sprintf("%.2f", filter(stats_res, par == "kappa", stat == "correlation")$median[1:4]),
  ) +
  theme_ggdist() +
  theme(
    panel.spacing = unit(1, "lines"),
    axis.text.x = element_text(
      angle = 80,
      hjust = 0.5,
      vjust = 0.5
    )
  )
)


(p2 <- stats %>%
  filter(par == "kappa", stat == "RMSE") %>%
  ggplot(aes(value, y = Nobs)) +
  geom_segment(
    data = filter(stats_res, par == "kappa", stat == "RMSE"),
    aes(x = median, yend = -Inf),
    color = "gray65",
    linetype = "dashed"
  ) +
  stat_halfeye(normalize = "groups", fill = hue_pal()(2)[1]) +
  facet_grid(cols = vars(par),
             scales = "free_x") +
  scale_x_continuous(
    "RMSE (logarithmic spacing)",
    breaks = filter(stats_res, par == "kappa", stat == "RMSE")$median,
    labels = sprintf("%.2f", filter(stats_res, par == "kappa", stat == "RMSE")$median),
    transform = "log"
  ) +
  theme_ggdist() +
  theme(
    panel.spacing = unit(1, "lines"),
    axis.text.x = element_text(
      angle = 80,
      hjust = 0.5,
      vjust = 0.5
    )
  )
)


(p3 <- stats %>%
  filter(par == "pmem", stat == "correlation") %>%
  ggplot(aes(value, y = Nobs)) +
  geom_segment(
    data = filter(stats_res, par == "pmem", stat == "correlation"),
    aes(x = median, yend = -Inf),
    color = "gray65",
    linetype = "dashed"
  ) +
  stat_halfeye(normalize = "groups", fill = hue_pal()(2)[2]) +
  facet_grid(cols = vars(par),
             scales = "free_x") +
  scale_x_continuous(
    "Correlation (regular spacing)",
    breaks = filter(stats_res, par == "pmem", stat == "correlation")$median,
    labels = sprintf("%.2f", filter(stats_res, par == "pmem", stat == "correlation")$median)
  ) +
  theme_ggdist() +
  theme(
    panel.spacing = unit(1, "lines"),
    axis.text.x = element_text(
      angle = 80,
      hjust = 0.5,
      vjust = 0.5
    )
  )
)

(p4 <- stats %>%
  filter(par == "pmem", stat == "RMSE") %>%
  ggplot(aes(value, y = Nobs)) +
  geom_segment(
    data = filter(stats_res, par == "pmem", stat == "RMSE"),
    aes(x = median, yend = -Inf),
    color = "gray65",
    linetype = "dashed"
  ) +
  stat_halfeye(normalize = "groups", fill = hue_pal()(2)[2]) +
  facet_grid(cols = vars(par),
             scales = "free_x") +
  scale_x_continuous(
    "RMSE (logarithmic spacing)",
    breaks = filter(stats_res, par == "pmem", stat == "RMSE")$median,
    labels = sprintf("%.2f", filter(stats_res, par == "pmem", stat == "RMSE")$median),
    transform = "log"
  ) +
  theme_ggdist() +
  theme(
    panel.spacing = unit(1, "lines"),
    axis.text.x = element_text(
      angle = 80,
      hjust = 0.5,
      vjust = 0.5
    )
  )
)

(p1 + p2) / (p3 + p4)










