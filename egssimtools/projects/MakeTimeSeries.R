rm(list = ls())

# Show speciation metrics through time

library(egssimtools)
library(tidyverse)
library(patchwork)
library(ggsim)

# Read the data
root <- "/media/raphael/bigass/simulations/EGS/mutator/"
data <- readRDS(paste0(root, "simulations.rds"))

head(data)

# Plot one variable through time across h and s,
# for one value of sigmaI
plot_this <- function(data, variable, ylab, title, color) {

  p <- data %>%
    filter(ecosel <= 3) %>%
    mutate(time = time / 1000) %>%
    mutate(ecosel = fct_rev(factor(ecosel))) %>%
    group_by(sim) %>%
    mutate(last = last(get(variable))) %>%
    ungroup() %>%
    gglineplot(x = "time", y = variable, line = "sim") +
    xlab("Time (\U02715 1,000 generations)") +
    ylab(ylab) +
    ggtitle(title) +
    aes(color = last) +
    scale_color_gradient(low = "black", high = color) +
    labs(color = variable)
  p %>% facettize(rows = "ecosel", cols = "hsymmetry", prepend = c("s = ", "h = "))

}

# Combine plots from multiple variables to make a figure, for a given value of the
# sigmaI parameter
combine_plots <- function(data, sigmaI) {

  figname <- sprintf("timeseries_sigmaI_%s.png", sigmaI)

  variables <- c("EI", "RI", "SI")
  ylabs <- c("Ecological divergence", "Reproductive isolation", "Spatial isolation")
  titles <- ylabs %>% map2_chr(c("A)", "B)", "C)"), ~ paste(.y, .x))
  colors <- c("lightgreen", "lightblue", "coral")

  args <- list(list(data), variables, ylabs, titles, colors)
  p <- args %>% pmap(plot_this)
  fig <- p[[1]] | p[[2]] | p[[3]]
  ggsave(figname, fig, width = 15, height = 9, dpi = 400)
  return (fig)

}

# Make and save figures for each value of the meta-parameter
data %>%
  group_by(scaleI1) %>%
  nest() %>%
  mutate(figure = map2(data, scaleI1, combine_plots))
