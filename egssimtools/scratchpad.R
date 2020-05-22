rm(list = ls())

library(egssimtools)
library(tidyverse)
library(tidygraph)
library(ggraph)

root <- "/media/raphael/bigass/simulations/EGS/genomes"
roots <- fetch_dirs(root, level = 2, pattern = "sim")
arch <- read_genome_architecture(roots[[1]])

arch <- read_network_architecture(roots[[1]])
arch

arch

col_vls <- c("forestgreen", "goldenrod", "grey")
col_lbs <- 0:2
col_nm <- "Trait"

ggraph(arch, layout = "graphopt") +
  geom_edge_link(aes(color = factor(trait), alpha = abs(weight))) +
  geom_node_point(aes(fill = factor(trait), size = degree), shape = 21) +
  scale_edge_color_manual(values = col_vls, labels = col_lbs, name = col_nm) +
  scale_fill_manual(values = col_vls, labels = col_lbs, name = col_nm) +
  labs(size = "Degree", edge_alpha = "|Weight|") +
  theme_void()


