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

data <- read_data(
  root, c("time", "genome_Fst"), by = c(1, 1), dupl = c(300, 1),
  architecture = TRUE, archfile = archfile
)

data <- data %>% mutate(locus = factor(locus))
nloci <- length(unique(data$locus))



# Plot locus Fst through time for each trait
p <- gglineplot(
  data, x = "time", y = "genome_Fst", line = "locus",
  mapping = aes(color = trait)
)
p <- facettize(p, rows = "trait", prepend = "trait ")
p
p + aes(color = trait)
