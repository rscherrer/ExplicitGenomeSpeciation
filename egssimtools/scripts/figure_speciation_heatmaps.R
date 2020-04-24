rm(list = ls())

library(egssimtools)
library(tidyverse)
library(cowplot)

# Settings
variables <- c("EI", "RI", "SI")
colors <- c("lightgreen", "lightblue", "coral")

# Read the data
data <- readRDS("data/population_wide_data.rds")

# Little extra tweaks
data$scaleI <- factor(data$scaleI, levels(data$scaleI)[c(1, 3, 2)])
levels(data$scaleI) <- c("0%", "50%", "100%")
data$mutation <- factor(data$mutation)
levels(data$mutation) <- c("0.0001", "0.001")

smr <- data %>%
  mutate(mutation = fct_rev(factor(mutation))) %>%
  group_by(hsymmetry, ecosel, dispersal, mutation, scaleA, scaleI, simulation) %>%
  summarize(EI = last(EI), SI = last(SI), RI = last(RI)) %>%
  group_by(hsymmetry, ecosel, dispersal, mutation, scaleA, scaleI) %>%
  summarize(EI = mean(EI), SI = mean(SI), RI = mean(RI)) %>%
  ungroup()

# Our heatmap figures
p <- mapply(function(variable, color) {

    ggplot(smr, aes(x = hsymmetry, y = ecosel, fill = get(variable))) +
      facet_grid(
        mutation ~ scaleI,
        labeller = do.call("labeller", make_facet_labels(
          data, c("mutation", "scaleI"),
          prefixes = c(mutation = "\U003BC", scaleI = "\U003C3\U0026A")
        ))) +
      geom_tile() +
      theme_bw() +
      ylab("Divergent selection") +
      xlab("Habitat symmetry") +
      labs(fill = variable) +
      scale_fill_gradient(low = "black", high = color)

}, variables, colors, SIMPLIFY = FALSE)

p[[1]]
p[[2]]
p[[3]]

fig <- plot_grid(plotlist = p, nrow = 3, labels = c("A", "B", "C"))
fig
ggsave("figures/heatmaps.png", width = 4, height = 7, dpi = 300)
