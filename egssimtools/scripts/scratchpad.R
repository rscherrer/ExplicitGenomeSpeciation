rm(list = ls())

library(tidyverse)
library(egssimtools)
library(patchwork)

root <- "data"

roots <- fetch_dirs(root, "example", level = 1)

read_data(roots[1], c("time", "EI"))
read_sim(roots[1], "EI")
read_pop(roots[1], "individual_trait", by = 3)
read_genome(roots[1], "genome_Fst", architecture = TRUE)
read_network(roots[1], "network_corbreed", architecture = TRUE)
