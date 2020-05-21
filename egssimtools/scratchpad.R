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
