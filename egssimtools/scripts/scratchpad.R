rm(list = ls())

library(tidyverse)
library(egssimtools)
library(patchwork)

root <- "data/example_1"

nedges <- guess_nedges(root)
read_data(root, c("time", "network_corbreed"), dupl = c(nedges, 1))

read_arch(root)
read_arch_genome(root)
read_arch_network(root, as_list = FALSE) %>% as.list %>% str
