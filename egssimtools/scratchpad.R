rm(list = ls())

library(egssimtools)
library(tidyverse)
library(ggsim)

root <- "/media/raphael/bigass/simulations/EGS/genomes3/"
roots <- fetch_dirs(root, "sim", level = 2)

root <- roots[150]
archfile <- root %>%
  str_replace("^.*arch", "arch") %>%
  str_replace("txt.*$", "txt")



data <- read_loci(root, "genome_Fst", architecture = TRUE, archfile = archfile)
data <- data %>% mutate(locus = factor(locus))
nloci <- length(unique(data$locus))

# Plot densities through time


# Plot lines
data <- smoothen_data(data, x = "time", y = "genome_Fst", line = "locus", span = 0.2)
p <- gglineplot(
  data, x = "time", y = "genome_Fst", line = "locus",
  mapping = aes(color = trait)
)
p <- facettize(p, rows = "trait", prepend = "trait ")
p
