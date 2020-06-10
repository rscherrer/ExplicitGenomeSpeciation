# In this script we test the functions of the package

rm(list = ls())

library(egssimtools)
library(tidyverse)
library(ggsim)

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
plot_trait_histogram(simulation, trait = 1, t = 2900, color = "seagreen", is_density = TRUE)
plot_genome_scan(simulation, "genome_Fst", t = 2900, colvar = "trait", col_as_factor = TRUE)
plot_genome_distribution(simulation, "genome_Fst", graph = list(grouping = "trait", plot_type = "boxplot"))
plot_genome_heatmap(simulation, "genome_Fst")

# Higher-level functions

root <- "/media/raphael/bigass/simulations/EGS/EGS_sim1"

find_extinct(root, level = 1)
find_missing(root, level = 1)
find_extant(root, level = 1)
find_completed(root, level = 1)
collect_status(root, level = 1)
collect_parameters(root, level = 1, parnames = c("ecosel", "hsymmetry"))
collect_simulations(root, level = 1, "RI", parnames = c("ecosel", "hsymmetry"))

