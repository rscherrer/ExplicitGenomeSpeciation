# In this script we test the functions of the package

rm(list = ls())

library(egssimtools)
library(tidyverse)

###

# Can I launch a simulation from R?

#
# Need to be in the target directory
#

cd

###

simulation <- "../build/release"

# Low-level functions for single mutations

read_paramfile(paste0(simulation, "/paramlog.txt"))
read_parameters(simulation)
read_binary(paste0(simulation, "/time.dat"))
read_data(simulation, "time")
read_time(simulation)
read_population_size(simulation)
read_population(simulation, "EI")
is_extinct(simulation)
is_missing(simulation)
read_archfile(paste0(simulation, "/architecture.txt"))
read_architecture(simulation)
read_individuals(simulation, "individual_trait")
read_ecotypes(simulation)
read_habitats(simulation)
read_loci(simulation, "genome_Fst")
read_edges(simulation, "network_corgen")
read_resources(simulation)
read_means(simulation)
read_ecotype_means(simulation)
read_ecotype_sizes(simulation)
plot_trait_density(simulation, trait = 1)
plot_trait_density2D(simulation, traits = c(2, 1), t = 100)


# Higher-level functions

root <- "/media/raphael/bigass/simulations/EGS/EGS_sim1"

find_extinct(root, pattern = "sim_")
find_missing(root, pattern = "sim_")
simulations <- find_extant(root, pattern = "sim_", pb = FALSE)
collect_parameters(simulations[1:4], parnames = c("ecosel", "hsymmetry"))
collect_simulations(simulations[1:4], "RI", parnames = c("ecosel", "hsymmetry"))

