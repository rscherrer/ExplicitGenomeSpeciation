# Sampling of parameter space for looking into genomic divergence

rm(list = ls())

library(egssimtools)
library(tidyverse)
library(ggsim)
library(ggnewscale)
library(cowplott)

data <- readRDS("data/simulations.rds")

set.seed(3)

data <- data %>%
  filter(mutation == 0.0001 & scaleI == "0 0 0") %>%
  shrink(c("EI", "RI", "SI"), list(c("simulation", "ecosel", "hsymmetry")), how = last) %>%
  mutate(completed = EI > 0.9) %>%
  group_by(hsymmetry, completed) %>%
  mutate(picked = seq(n()) %in% sample(n(), 5))

p1 <- data %>% ggplot(aes(x = RI, y = SI, color = EI)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  coord_fixed(xlim = c(-0.1, 1), ylim = c(-0.1, 1)) +
  scale_color_gradient(low = "black", high = "lightgreen") +
  guides(size = FALSE, alpha = FALSE) +
  labs(
    x = "Reproductive isolation",
    y = "Spatial isolation",
    color = "Ecological divergence",
    fill = "Habitat symmetry"
  ) +
  geom_point(
    data = data %>% filter(picked),
    aes(fill = hsymmetry),
    color = "black", shape = 21, size = 2.5
  ) +
  scale_fill_gradient(low = "coral", high = "brown4")

p2 <- data %>%
  ggheatmap("EI", "hsymmetry", "ecosel", how = mean) +
  scale_fill_gradient(low = "black", high = "lightgreen") +
  new_scale("fill") +
  geom_point(data = data %>% filter(picked), aes(fill = hsymmetry), shape = 21, size = 2.5) +
  scale_fill_gradient(low = "coral", high = "brown4") +
  labs(x = "Habitat symmetry", y = "Divergent selection")

plot_grid(
  p1 + theme(legend.position = "none"),
  p2 + theme(legend.position = "none"),
  labels = c("A", "B")
) %>%
  plot_grid(
    get_legend(p1 + theme(legend.position = "top")),
    ., nrow = 2, rel_heights = c(1, 4)
  )
ggsave("figures/sampling.png", width = 6, height = 4, dpi = 300)

data %>%
  ungroup() %>%
  filter(picked, !completed) %>%
  select(simulation) %>%
  droplevels() %>%
  unlist() %>%
  first()
