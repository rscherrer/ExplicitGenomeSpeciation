rm(list = ls())

library(egssimtools)
library(tidyverse)

root <- "/media/raphael/bigass/simulations/EGS/EGS_sim1"

roots <- fetch_dirs(root, pattern = "sim_", level = 1)[1:10]

data <- read_data(roots[1], c("time", "EI", "RI", "SI"))
pars <- read_parameters(roots[1], combine = FALSE, flatten = TRUE, as_numeric = "scaleA")
pars$scaleA1
