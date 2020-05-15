rm(list = ls())

library(egssimtools)
library(tidyverse)
library(patchwork)

# Load the data from the experiment
figname <- "heatmaps.png"
root <- "/media/raphael/bigass/simulations/EGS/mutator/"
data <- readRDS(paste0(root, "simulations.rds"))

head(data)

# Heatmap of the final value of the variables across parameter space
plot_this <- function(variable, color, title) {
  p <- ggheatmap(
    data, variable = variable, x = "hsymmetry", y = "ecosel", reduce = "sim",
    how = c(last, mean), keep = c("scaleI1")
  ) +
    scale_fill_gradient(low = "black", high = color) +
    xlab("Habitat symmetry") +
    ylab("Ecological trade-off") +
    ggtitle(title)
  p %>% facettize(cols = "scaleI1", prepend = "sigma[I]==", parsed = "scaleI1")
}

# Assemble the plots for each variable into one figure
variables <- c("EI", "RI", "SI")
colors <- c("lightgreen", "lightblue", "coral")
titles <- c("A) Ecological divergence", "B) Reproductive isolation", "C) Spatial isolation")

args <- list(variables, colors, titles)

p <- args %>% pmap(plot_this)
fig <- p[[1]] / p[[2]] / p[[3]]

ggsave(figname, fig, width = 6, height = 8, dpi = 300)
