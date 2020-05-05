# Use this script to make a database summarizing data over many simulations

rm(list = ls())

library(egssimtools)
library(tidyverse)

# Make a data base of simulations through time

root <- "/media/raphael/bigass/simulations/EGS/mutator"
variables <- c("time", "EI", "RI", "SI")
parnames <- c("ecosel", "hsymmetry", "mutation", "scaleI")
as_numeric <- c("ecosel", "hsymmetry", "mutation", "scaleI")

roots <- fetch_dirs(root, level = 2, pattern = "sim_")
roots <- find_extant(roots)
data <- roots %>% map_dfr(read_data, variables, parnames = parnames, as_numeric = as_numeric, .id = "simulation")

head(data)

saveRDS(data, paste0(root, "/simulations.rds"))
