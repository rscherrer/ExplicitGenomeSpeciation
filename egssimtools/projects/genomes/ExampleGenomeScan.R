rm(list = ls())

library(egssimtools)
library(tidyverse)
library(ggsim)
library(patchwork)

# This is saving a hundred genomes

root <- "/media/raphael/bigass/simulations/EGS/genomes/"
data <- readRDS(paste0(root, "simulations.rds"))
backup <- data
data <- backup

head(data)

# move to mkdatabase
data <- data %>% mutate(locus = rep(1:300, nrow(data) / 300))
data$trait <- factor(data$trait)
data <- data %>% filter(time == 19900)

# Plot a single replicate of a single architecture across parameter space
# for the last time point only

# Function to plot a genome scan for a given simulation
plot_this <- function(data) {

  colors <- c("forestgreen", "goldenrod", "grey")

  # Plot Fst along the genome
  p <- data %>%
    ggplot(aes(x = location, y = genome_Fst, color = factor(trait))) +
    geom_point(alpha = 0.8) +
    theme_bw() +
    scale_color_manual(values = colors) +
    labs(x = "Location", y = "Fst", color = "Trait")

  # Plot the distribution of Fst across traits
  p2 <- data %>%
    ggplot(aes(x = trait, y = genome_Fst, color = trait)) +
    geom_violin() +
    theme_bw() +
    scale_color_manual(values = colors) +
    labs(x = "Location", y = "Fst", color = "Trait")

  # Combine both into one example figure for that specific architecture
  fig <- p + theme(legend.position = "none") | p2
  fig <- fig +
    plot_layout(widths = c(3, 1)) +
    plot_annotation(tag_levels = "A")

  fig

}

# Make a plot per simulation
newdata <- data %>%
  group_by(archfile, seed, hsymmetry, scaleI1) %>%
  nest() %>%
  mutate(figure = map(data, plot_this))

newdata$figure[[2]]

# Saving function
save_this <- function(figure, figname) {

  ggsave(figname, figure, width = 9, height = 3, dpi = 300)

}

# Save the figures
figname <- "genomescan_%s_h_%s_sigmaI_%s_seed_%s.png"
newdata %>%
  mutate(arch = gsub(".txt", "", archfile)) %>%
  mutate(figname = sprintf(figname, arch, hsymmetry, scaleI1, seed)) %>%
  mutate(saved = walk2(figure, figname, save_this))




