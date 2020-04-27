# In this script we test the functions of the package

rm(list = ls())

library(egssimtools)
library(tidyverse)

###

# We see that when the mutation rate is high, reproductive isolation in asymmetric
# habitats evolves even less than when mutation is high. Why is that? To answer
# we compare the trait distributions through time in high and low mutation scenarios...
# Hypothesis: the high mutation rate increases the variance in ecological trait
# such that ecotypes are not ecologically isolated enough for RI to be advantageous.
# If this is true, we should see RI evolving at higher mutation rates, where ecotypes
# become more isolated.

root <- "/media/raphael/bigass/simulations/EGS_extra/EGS_sim3"

data <- collect_simulations(
  root, variables = c("EI", "RI", "SI"), parnames = c("ecosel", "hsymmetry"),
  to_numeric = c("ecosel", "hsymmetry")
)

# Can I launch a simulation from R?

#
# Need to be in the target directory
#

data <- readRDS("data/population_wide_data.rds")
head(data)

where <- "../cluster/test"
exe <- "../EGS"
parfile <- NULL

lines <- c(
  "#!/bin/bash",
  paste("cd", where),

)

jobfile <- file("job.sh")
writeLines(c("#!/bin/bash", "cd ../cluster/test", '../EGS'), jobfile)
close(jobfile)

system("chmod u+x ./job.sh")
system("./job.sh")
system("ls ../cluster/test")

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

