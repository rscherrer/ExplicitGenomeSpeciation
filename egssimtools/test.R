# In this script we test the functions of the package

rm(list = ls())

simulation <- "../build/release"

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
find_extinct("../build", pattern = "release")
find_missing("../build", pattern = "release")

read_means <- function(folder) {

  data <- read_data(folder, "resources")
  t <- read_time(folder)
  t <- rep(t, each = 4) # per time point
  lapply(split(data, t), split, c(1, 1, 2, 2)) # per habitat

}

read_resources(simulation)
read_means(simulation)


