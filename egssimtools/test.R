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

read_individuals(simulation)
