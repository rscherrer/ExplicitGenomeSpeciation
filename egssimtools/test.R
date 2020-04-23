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
simulations <- find_extant(root, pattern = "sim_", pb = FALSE)
collect_parameters(simulations, parnames = c("ecosel", "hsymmetry"))
collect_simulations(root, "RI", parnames = c("ecosel", "hsymmetry"), filters = "ecosel == 1")

# High-level functions: produce figures and overall results over multiple simulations

root <- "/media/raphael/bigass/simulations/EGS/EGS_sim1"

add_summaries <- function(data) data %>% group_by(simulation) %>% mutate(x = last(RI)) %>% ungroup() %>% select(x)

plot_simulations(root, "RI", facet_rows = "hsymmetry")

variable <- "RI"
facet_rows <- NULL
facet_cols <- "hsymmetry"
facet_wrapped <- FALSE
label_facets <- FALSE
reverse_order <- NULL
facet_prefixes <- NULL
color_by <- NULL
color_by_numeric <- TRUE
colors <- NULL
pb <- TRUE
verbose <- TRUE
sep <- " "
add_summaries <- NULL
filters <- NULL
xvariable <- NULL
pattern = "sim_"

# broken
plot_simulations(
  root,
  "RI",
  facet_cols = "hsymmetry",
  facet_rows = "ecosel",
  reverse_order = "ecosel",
  color_by = "x",
  xvariable = "EI",
  add_summaries = function(data) data %>% group_by(simulation) %>% mutate(x = last(EI)) %>% ungroup() %>% select(x)
)

simulations <- find_extant(root)

data <- collect_simulations(simulations[1:4], c("RI", "EI", "SI"), parnames = c("ecosel", "hsymmetry"))
head(data)
