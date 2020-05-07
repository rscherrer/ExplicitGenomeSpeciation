# Read the simulations and assemble them into a data set

rm(list = ls())

library(egssimtools)
library(tidyverse)
library(ggsim)

root <- "/media/raphael/bigass/simulations/EGS/genomes/"

variables <- c("time", "Fst", "Gst", "Qst", "Cst", "varG", "varA", "varN", "varT")
by <- c(1, 3, 3, 3, 3, 3, 3, 3, 3)
dupl <- c(1, 1, 1, 1, 1, 1, 1, 1, 1)
parnames <- c("hsymmetry", "ecosel", "scaleI", "seed")
as_numeric <- parnames

data <- collect_sims(
  root, variables, by = by, dupl = dupl, parnames = parnames,
  as_numeric = as_numeric, check_extant = FALSE, level = 1, pattern = "sim"
)

data %>%
  gglineplot(x = "time", y = "Cst3", line = "sim") %>%
  ggfacet("hsymmetry", "scaleI1", prepend = c(hsymmetry = "h = ", scaleI1 = "sigma_i = "))

data <- collect_sims(root, c("time", "genome_Cst"), by = c(1, 1), dupl = c(300, 1),
                     parnames = parnames, check_extant = FALSE, level = 1, pattern = "sim")

head(data)


