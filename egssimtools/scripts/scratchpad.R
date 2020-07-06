rm(list = ls())

library(egssimtools)

root <- "/media/raphael/bigass/simulations/EGS/genomes6/batch4/sim_ecosel_1.4_hsymmetry_1_scaleA_0_0_0_scaleI_1_1_1_archfile_arch1.txt_r2/"

filename <- file.path(root, "freezer.dat")
genome <- read_bitset(filename)

nloci <- guess_nloci(root)

# Number of zeros added at the end of each individual genome
nextra <- 64 - ((2 * nloci) %% 64)

# Number of bits per individual genome
nbits <- (2 * nloci) + nextra

# Number of individuals
ninds <- length(genome) / nbits

# Make a large dataset
inds <- rep(seq(ninds), each = nbits)
loci <- rep(c(rep(seq(nloci), 2), rep(NA, nextra)), ninds)
haps <- rep(c(rep(seq(2), each = nloci), rep(NA, nextra)), ninds)
data <- tibble::tibble(
  allele = genome,
  locus = loci,
  ind = inds,
  hap = haps
)

# Remove the extra zeros (they were added to convert individual genomes
# into 64bit integers)
data <- data %>% drop_na()

# Convert to the wide format (don't thinks we will need the long-format
# in this study) = locus-wise
data <- data %>% mutate(hap = stringr::str_replace(hap, "^", "hap"))
data <- data %>% pivot_wider(names_from = "hap", values_from = "allele")

t <- 19900

# Read the time point when each individual lived
popsizes <- read_sim(root, "population_size")
popsizes <- popsizes %>% dplyr::filter(time %in% t)
times <- do.call("c", with(popsizes, map2(time, population_size, ~ rep(.x, .y))))
data$time <- rep(times, each = nloci)

head(data)

# Derive additional data from the alleles
data <- data %>% mutate(allcount = hap1 + hap2)

# Add additional data if needed
extra <- c("individual_ecotype")

extradata <- read_pop(root, extra)
