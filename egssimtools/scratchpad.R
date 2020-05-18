rm(list = ls())

library(egssimtools)
library(tidyverse)
library(ggsim)
library(patchwork)

root <- "data"

# Fetch simulations folders
fetch_dirs(roots = root, level = 1)

root <- "data/example_1"
variables <- c("time", "genome_Fst")
data <- read_data(root, variables, dupl = c(300, 1), architecture = TRUE)
head(data)

# Genome scan at a certain time
t <- 19900
data <- data %>% filter(time == t)
ggplot(data, aes(x = location, y = genome_Fst)) +
  geom_point() +
  theme_bw()

dplot_genome_scan(root, y = "genome_Fst")
dplot_genome_heatmap(root, "genome_Fst")
dplot_genome_violin(root, y = "genome_Fst", x = "trait")

dplot_population_density(root, y = "individual_trait", by = 3, j = 1)
dplot_population_bin2d(root, y = "individual_trait", by = 3, j = 1)

dplot_simulation_line(root, y = "EI", by = 1, j = 1)
