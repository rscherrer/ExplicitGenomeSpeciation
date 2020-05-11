rm(list = ls())

library(egssimtools)
library(tidyverse)

root <- "/media/raphael/bigass/simulations/EGS/mutator/"

variables <- c("time", "EI", "RI", "SI")
parnames <- c("hsymmetry", "ecosel", "scaleA", "scaleI")

data <- collect_sims(root, level = 2, pattern = "sim_", variables = variables,
                     parnames = parnames, check_extant = TRUE,
                     as_numeric = parnames)
data
