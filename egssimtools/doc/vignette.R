## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message = FALSE---------------------------------------------------------
library(egssimtools)
library(tidyverse)
library(cowplot)

## -----------------------------------------------------------------------------
root <- "../data/example_1"
data <- read_data(
  root, variables = c("time", "EI", "RI", "SI"), 
  parnames = c("ecosel", "hsymmetry")
)
head(data, 4)

## -----------------------------------------------------------------------------
data <- read_data(root, variables = c("time", "Fst"), by = c(1, 3))
head(data, 4)

## -----------------------------------------------------------------------------
data <- read_data(root, variables = c("time", "genome_Fst"), by = c(1, 300))
head(data[, 1:6], 4)

## -----------------------------------------------------------------------------
data <- read_data(
  root, variables = c("time", "genome_Fst", "genome_freq", "genome_varA"), 
  dupl = c(300, 1, 1, 1)
)
head(data, 4)

## -----------------------------------------------------------------------------
data <- read_data(
  root, variables = c("time", "individual_trait"), by = c(1, 3),
  dupl = list("population_size", 1)
)
head(data, 4)

## -----------------------------------------------------------------------------
read_parameters(root, c("ecosel", "hsymmetry"))

## ---- message = FALSE---------------------------------------------------------
data <- collect_sims(
  root = "../data", 
  variables = c("time", "EI"), 
  parnames = c("ecosel",  "hsymmetry"), 
  id_column = "sim", 
  level = 1,
  pattern = "example",
  verbose = FALSE)
head(data, 4)

## -----------------------------------------------------------------------------
arch <- read_architecture(root)
str(arch)

## -----------------------------------------------------------------------------
loci <- read_genome_architecture(root)
head(loci, 4)

## ---- fig.width = 4-----------------------------------------------------------
dplot_genome_scan(root, y = "genome_Fst")

## ---- fig.width = 5-----------------------------------------------------------
dplot_genome_scan(root, y = "genome_Fst", t = 19900) + 
  aes(color = factor(trait)) +
  labs(x = "Locus", y = parse(text = "F[ST]"), color = "Trait") +
  scale_color_manual(values = c("forestgreen", "goldenrod", "lightgrey"))

## ---- fig.width = 5-----------------------------------------------------------
dplot_genome_heatmap(root, y = "genome_Fst") +
  scale_fill_gradient(low = "darkred", high = "gold") +
  labs(x = "Time (generations)", y = "Locus", fill = parse(text = "F[ST]"))

## ---- fig.width = 4-----------------------------------------------------------
dplot_genome_violin(root, y = "genome_Fst", x = "trait") +
  aes(color = factor(trait)) +
  scale_color_manual(values = c("forestgreen", "goldenrod", "lightgrey")) +
  labs(x = "Trait", y = parse(text = "F[ST]"), color = "Trait")

## ---- fig.width = 3-----------------------------------------------------------
dplot_population_density(
  root, y = "individual_trait", by = 3, j = 1, fill = "lightgreen"
) +
  labs(x = "Ecological trait", y = "Density") +
  xlim(c(-2, 2))

## ---- fig.width = 5-----------------------------------------------------------
dplot_population_bin2d(
  root, y = "individual_trait", by = 3, j = 1, bins = 100
) +
  labs(x = "Time (generations)", y = "Ecological trait", fill = "Count") +
  scale_fill_continuous(type = "viridis")

## ---- fig.width = 6-----------------------------------------------------------
p1 <- dplot_simulation_line(root, y = "Fst", by = 3, j = 1) +
  labs(x = "Time (generations)", y = parse(text = "F[ST]"))
p2 <- dplot_simulation_line(root, y = "RI", x = "EI") +
  labs(x = "Ecological divergence", y = "Reproductive isolation")
plot_grid(p1, p2, labels = c("A", "B"))

