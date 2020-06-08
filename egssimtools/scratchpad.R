rm(list = ls())

library(egssimtools)

root <- system.file("extdata", "example_1", package = "egssimtools")
plot_genome_density(root, "genome_Fst", x = "trait", mapping = list(fill = "trait"))

