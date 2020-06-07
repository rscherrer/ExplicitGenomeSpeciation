rm(list = ls())

library(egssimtools)

root <- "/media/raphael/bigass/simulations/EGS/genomes3/"
roots <- fetch_dirs(root, pattern = "sim_", level = 2)
root <- roots[1]

read_data(root, c("time", "Fst"), by = c(1, 3))
read_sim(root, "Fst", by = 3)
read_indiv(root, "individual_trait")

read_edges()
