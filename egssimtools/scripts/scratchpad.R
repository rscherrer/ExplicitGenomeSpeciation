rm(list = ls())

library(egssimtools)

roots <- fetch_dirs("data", "example", 1)
find_extant(roots)
