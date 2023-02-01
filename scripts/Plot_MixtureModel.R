# start fresh
rm(list = ls()) # clean up work space
graphics.off()  # switch off graphics device

# load required packages
library(ggplot2)
library(here)

# load function to clean up plots
source(here("functions","clean_plot.R"))


# settings for mixture model plot
mus <- c(0,1.7,0)
kappas <- c(10^-3,7,7)
weights <- c(0.5,0.2,1)
colors  <- c("firebrick", "chartreuse4", "dodgerblue4","black")

# weighted density function for the vonMises distribution
weighted_vonMises_density <- function(x, mu, kappa, log = F, weight = 1){
  weighted_vonMises_density <- dvon_mises(x = x,mu = mu , kappa = kappa, log = log) * weight
  return(weighted_vonMises_density)
}

# sum of three von Mises densities for plotting the mixture density
sum_vonMises_density <- function(x,mus,kappas,log = F, weights) {
  sum_density <- dvon_mises(x, mu = mus[1], kappa = kappas[1])*(weights[1]/sum(weights)) +
    dvon_mises(x, mu = mus[2], kappa = kappas[2])*(weights[2]/sum(weights)) +
    dvon_mises(x, mu = mus[3], kappa = kappas[2])*(weights[3]/sum(weights))
  return(sum_density)
}

# define defaults for clean plots
clean_plot <- theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    axis.line.x = element_line(color = 'black'),
                    axis.line.y = element_line(color = 'black'),
                    axis.text.y = element_blank(),
                    legend.key = element_rect(fill = 'white'),
                    text = element_text(size = 15),
                    line = element_line(size = 1),
                    axis.ticks = element_line(size = 1))

# ggplot to illustrate the mixture model
plot_mixModel <- ggplot() +
  stat_function(fun = "weighted_vonMises_density",
                args = list(mu = mus[1], kappa = kappas[1], log = F, weight = weights[1]), 
                col = colors[1], lwd = 1.5, lty = "dashed") +
  stat_function(fun = "weighted_vonMises_density",
                args = list(mu = mus[2], kappa = kappas[2], log = F, weight = weights[2]),
                col = colors[2], lwd = 1.5, lty = "longdash") +
  stat_function(fun = "weighted_vonMises_density",
                args = list(mu = mus[3], kappa = kappas[3], log = F, weight = weights[3]), 
                col = colors[3], lwd = 1.5, lty = "solid") +
  stat_function(fun = "sum_vonMises_density",
                args = list(mus = mus, kappas = kappas, weights = weights), 
                col = colors[4], lwd = 1.5, lty = "dotted") +
  scale_x_continuous(limits = c(-pi,pi),
                     breaks = c(-3,-2,-1,0,1,2,3)) +
  annotate(
    "text",label = "guessing", x = -2, y = 0.15,
    colour = colors[1], size = 6) +
  annotate(
    "text",label = "other object", x = mus[2], y = 0.3,
    colour = colors[2], size = 6) +
  annotate(
    "text",label = "cued object", x = 0.9, y = 0.8,
    colour = colors[3], size = 6) +
  annotate(
    "text", label = "Mixture", x = 0, y = 0.3,
    colour = colors[4], size = 6) +
  labs(x = "deviation from cued object (in radians)",
       y = "probability of report") +
  # +
  clean_plot

# show plot
plot_mixModel

# save plot with high resolution
ggsave(
  filename = here("figures/plot_mixModel.jpeg"),
  plot = plot_mixModel, width = 8, height = 6
)
