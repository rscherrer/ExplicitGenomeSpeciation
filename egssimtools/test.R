# In this script we test the functions of the package

rm(list = ls())

library(egssimtools)
library(tidyverse)

simulation <- "../build/release"

# Low-level functions: work on single simulations

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


# Medium-level functions: can combine and process data from different simulations but the output needs to be processed to produce full-fledged figures and results

root <- "/media/raphael/bigass/simulations/EGS/EGS_sim1"

find_extinct(root, pattern = "sim_")
find_missing(root, pattern = "sim_")


# High-level functions: produce figures and overall results over multiple simulations

root <- "/media/raphael/bigass/simulations/EGS/EGS_sim1"

plot_simulations(root, "RI", color_by = "hsymmetry", colors = c("lightgrey", "darkblue", "black", "yellow", "red"), color_by_numeric = FALSE, facet_rows = "ecosel", facet_cols = "hsymmetry", label_facets = TRUE, facet_prefixes = c("s", "h"), reverse_order = c("ecosel", "hsymmetry"), verbose = TRUE, pb = FALSE)
plot_simulations(root, "EI", facet_rows = "hsymmetry", facet_wrapped = TRUE, color_by = "ecosel", label_facets = TRUE, facet_prefixes = "s", verbose = TRUE, pb = FALSE)
