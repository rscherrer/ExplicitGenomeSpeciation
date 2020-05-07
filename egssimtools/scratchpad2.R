# Read the simulations and assemble them into a data set

rm(list = ls())

library(egssimtools)
library(tidyverse)

root <- "/media/raphael/bigass/simulations/EGS/genomes/"

root <- list.dirs(root, recursive = FALSE)[1]
variables <- c("time", "Fst", "Gst", "Qst", "Cst")
by <- c(1, 3, 3, 3, 3)
dupl <- c(1, 1, 1, 1, 1)
parnames <- c("hsymmetry", "ecosel", "scaleI", "seed")
as_numeric <- parnames


data <- collect_sims(
  root, variables, by = by, dupl = dupl, parnames = parnames,
  as_numeric = as_numeric, check_extant = FALSE, level = 1, pattern = "sim"
)


