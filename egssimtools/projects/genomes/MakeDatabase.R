rm(list = ls())

# Make a database of locus-specific metrics for all loci across all time points

library(egssimtools)

root <- "/media/raphael/bigass/simulations/EGS/genomes/"

# Load the data

variables <- c("Fst", "Gst", "Qst", "Cst", "varG", "varA", "varN")
variables <- paste0("genome_", variables)
variables <- c("time", variables)

by <- rep(1, length(variables))
dupl <- rep(1, length(variables))
dupl[1] <- 300

parnames <- c("hsymmetry", "scaleI", "seed", "archfile")

data <- collect_sims(
  root, level = 2, variables = variables, parnames = parnames, by = by,
  dupl = dupl, check_extant = FALSE, pattern = "sim"
)


# Add architecture metadata
data <- data %>%
  group_by(archfile) %>%
  nest() %>%
  mutate(arch = list(read_genome_architecture(root, filename = archfile))) %>%
  mutate(ntimes = nrow(data[[1]]) / nrow(arch[[1]])) %>%
  mutate(arch = map2(arch, ntimes, ~ map_dfr(seq(.y), function(i) .x))) %>%
  select(-ntimes) %>%
  unnest(cols = c(data, arch)) %>%
  ungroup()

saveRDS(data, paste0(root, "simulations.rds"))
