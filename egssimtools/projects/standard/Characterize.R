# Use this script to generate a heatmap across parameter space

rm(list = ls())

library(tidyverse)
library(egssimtools)

root <- "/media/raphael/bigass/simulations/EGS/standard/"
data <- readRDS(paste0(root, "/simulations.rds"))

head(data)

data <- data %>% mutate(scaleI1 = scaleI1 %>% factor %>% str_replace("^", "sigma[I]=="))

variables <- c("EI", "RI", "SI")
colors <- c("lightgreen", "lightblue", "coral")

plots <- list(variables, colors) %>% pmap(function(variable, color) {
  data %>%
    ggheatmap(
      variable, x = "hsymmetry", y = "ecosel", reduce = "simulation",
      how = c(last, mean), keep = "scaleI1"
    ) +
    facet_grid(. ~ scaleI1, labeller = label_parsed) +
    scale_fill_gradient(low = "black", high = color) +
    labs(x = "Habitat symmetry", y = "Ecological trade-off")
})

plot_grid(plotlist = plots, nrow = 3, labels = c("A", "B", "C"))

# ggsave("projects/standard/heatmaps.png", width = 5, height = 6, dpi = 300)

data <- data %>%
  mutate(
    time = time / 1000,
    hsymmetry = hsymmetry %>% factor %>% str_replace("^", "h = "),
    ecosel = ecosel %>% factor %>% str_replace("^", "s = ") %>% fct_rev
  )

labels <- c("Ecological divergence", "Reproductive isolation", "Spatial isolation")

plots <- data %>%
  group_by(scaleI1) %>%
  group_map(function(df, ...) {
    list(variables, colors, labels) %>%
      pmap(function(variable, color, label) {
        df %>%
          group_by(hsymmetry, ecosel, simulation) %>%
          mutate(last = last(get(variable))) %>%
          ungroup() %>%
          gglineplot(x = "time", y = variable, line = "simulation") +
          aes(color = last) +
          scale_color_gradient(low = "black", high = color) +
          facet_grid(ecosel ~ hsymmetry) +
          labs(x = "Time (thousands of generations)", y = label, color = variable)
      }) %>%
      plot_grid(plotlist = ., ncol = 3, labels = c("A", "B", "C"))
  })
ggsave("projects/standard/lineplots_additive.png", plots[[1]], width = 16, height = 9, dpi = 400)
ggsave("projects/standard/lineplots_intermediate.png", plots[[2]], width = 16, height = 9, dpi = 400)
ggsave("projects/standard/lineplots_epistatic.png", plots[[3]], width = 16, height = 9, dpi = 400)

plots <- data %>%
  group_by(scaleI1) %>%
  group_map(function(df, ...) {
    df %>%
      mutate(ecosel = ecosel %>% fct_rev) %>%
      gglineplot(x = "EI", y = "SI", line = "simulation") +
      aes(color = hsymmetry) +
      scale_color_manual(values = colorRampPalette(c("coral", "black"))(5)) +
      facet_wrap(. ~ ecosel) +
      labs(x = "Ecological divergence", y = "Spatial isolation", color = "Habitat symmetry")
  })
ggsave("projects/standard/phaseplots_ecospatial_additive.png", plots[[1]], width = 6, height = 5, dpi = 300)
ggsave("projects/standard/phaseplots_ecospatial_intermediate.png", plots[[2]], width = 6, height = 5, dpi = 300)
ggsave("projects/standard/phaseplots_ecospatial_epistatic.png", plots[[3]], width = 6, height = 5, dpi = 300)

plots <- data %>%
  group_by(scaleI1) %>%
  group_map(function(df, ...) {
    df %>%
      group_by(hsymmetry, ecosel, simulation) %>%
      mutate(last = last(RI)) %>%
      ungroup() %>%
      gglineplot(x = "EI", y = "RI", line = "simulation") +
      aes(color = last) +
      scale_color_gradient2(low = "purple", mid = "black", high = "lightblue") +
      facet_grid(ecosel ~ hsymmetry) +
      labs(x = "Ecological divergence", y = "Reproductive isolation", color = "RI")
  })
ggsave("projects/standard/phaseplots_ecorepro_additive.png", plots[[1]], width = 6, height = 8.5, dpi = 400)
ggsave("projects/standard/phaseplots_ecorepro_intermediate.png", plots[[2]], width = 6, height = 8.5, dpi = 400)
ggsave("projects/standard/phaseplots_ecorepro_epistatic.png", plots[[3]], width = 6, height = 8.5, dpi = 400)
