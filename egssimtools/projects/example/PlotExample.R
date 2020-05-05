
rm(list = ls())

library(egssimtools)
library(ggsim)
library(tidyverse)
library(cowplot)

root <- "/media/raphael/bigass/simulations/EGS/example"

# Read the simulation data in and save a database for easy access

roots <- fetch_dirs(root, level = 1, pattern = "sim_")
parnames <- c("ecosel", "hsymmetry", "scaleI", "seed")

data <- roots %>%
  map_dfr(
    read_data, c("time", "individual_trait"), by = c(1, 3), dupl = list("population_size", 1),
    parnames = parnames, as_numeric = parnames, .id = "simulation"
  )

saveRDS(data, paste0(root, "/simulations.rds"))

data <- readRDS(paste0(root, "/simulations.rds"))

# Plot densities of trait values through time

# Convert to long format
trait_cols <- colnames(data)[grep("individual_trait", colnames(data))]
data <- data %>%
  gather_("trait", "phenotype", trait_cols) %>%
  mutate(trait = trait %>% str_replace("individual_trait", "") %>% as.numeric)

# Facet labels
data <- data %>%
  mutate(
    hsymmetry = hsymmetry %>% str_replace("^", "h = "),
    ecosel = ecosel %>% str_replace("^", "s = "),
    mode = paste(hsymmetry, ", ", ecosel),
    scaleI1 = scaleI1 %>% str_replace("^", "sigma[I]==")
  )

# For each seed, make a figure with one panel per trait
# Each panel shows the evolution of the trait through time across four
# different genetic / ecological scenarios (in facets)

trait_labs <- c("Ecological trait", "Mating preference", "Neutral trait")
plotname <- "evoplot_seed_%s.png"
data %>%
  split(f = .[["seed"]]) %>%
  map(
    ~ .x %>% split(f = .[["trait"]]) %>%
      map(
        ~ ggplot(.x) +
          geom_bin2d(aes(x = time, y = phenotype), bins = 50) +
          facet_grid(mode ~ scaleI1, labeller = labeller(mode = label_value, scaleI1 = label_parsed)) +
          theme_bw() +
          scale_fill_continuous(type = "viridis") +
          labs(x = "Time (generations)", y = "Trait value") +
          ggtitle(trait_labs[.x$trait[1]]) +
          theme(legend.position = "none")
      ) %>%
      plot_grid(plotlist = ., ncol = 3, labels = c("A", "B", "C")) %>%
      ggsave(sprintf(plotname, .x$seed[1]), ., width = 10, height = 3, dpi = 400)
  )

