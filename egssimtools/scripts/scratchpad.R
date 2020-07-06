rm(list = ls())

library(egssimtools)

root <- "/media/raphael/bigass/simulations/EGS/genomes6/batch4/sim_ecosel_1.4_hsymmetry_1_scaleA_0_0_0_scaleI_1_1_1_archfile_arch1.txt_r2/"

filename <- file.path(root, "freezer.dat")

# How many bytes should there be?
# Number of loci: 300
# Number of alleles: 600
# Number of integers per individual: 10
n <- 600 %/% 64 + 1
# Generation: 19900
# Population size at that generation: 3534
data <- read_sim(root, "population_size")
size <- data %>% dplyr::filter(time == 19900)
size <- size$population_size
# Number of integers for the whole population:
nints <- n * size
# Number of bytes for the whole population:
nbytes <- nints * 8
nbytes

# And what is the actual size of the freezer?
file.size(filename)
