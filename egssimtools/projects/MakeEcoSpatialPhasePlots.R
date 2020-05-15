rm(list = ls())

# Show spatial isolation against ecological divergence

library(egssimtools)
library(tidyverse)
library(patchwork)
library(ggsim)

# Read the data
root <- "/media/raphael/bigass/simulations/EGS/mutator/"
data <- readRDS(paste0(root, "simulations.rds"))

head(data)

# Plot one variable against another, facetted across combinations of parameters
# For one value of the meta-parameter sigmaI
plot_this <- function(data, sigmaI) {

  title <- sprintf('"Parameter value " ~ sigma[I]==%s', sigmaI)

  p <- data %>%
    filter(ecosel <= 3) %>%
    group_by(sim) %>%
    mutate(last = last(RI)) %>%
    ungroup() %>%
    gglineplot(x = "EI", y = "SI", line = "sim") +
    xlab("Ecological divergence") +
    ylab("Spatial isolation") +
    ggtitle(parse(text = title)) +
    aes(color = hsymmetry) +
    scale_color_gradient(low = "coral", high = "black") +
    labs(color = "RI") +
    theme(legend.position = "none")
  p %>% facettize(rows = "ecosel", wrap = TRUE, prepend = "s = ")

}

# Apply the function to all values of the meta-parameter
newdata <- data %>%
  group_by(scaleI1) %>%
  nest() %>%
  mutate(plot = map2(data, scaleI1, plot_this))

newdata <- newdata[newdata$scaleI1 %>% order, ]

# Assemble the plots into one figure
fig <- newdata$plot %>% wrap_plots(ncol = 3, nrow = 1)
ggsave("phaseplot_ecospatial.png", fig, width = 13, height = 5, dpi = 400)
