rm(list = ls())

library(egssimtools)

# Load the data from the experiment
root <- "/media/raphael/bigass/simulations/EGS/standard/"

variables <- c("time", "EI", "RI", "SI")
parnames <- c("hsymmetry", "ecosel", "scaleA", "scaleI")

# This may take a while
data <- collect_sims(
  root, level = 2, pattern = "sim_", variables = variables,
  parnames = parnames, check_extant = TRUE, as_numeric = parnames
)

head(data)

# Save the data in a compressed format so we do not have to re-read them the next time
saveRDS(data, paste0(root, "simulations.rds"))
