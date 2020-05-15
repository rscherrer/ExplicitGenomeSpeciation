rm(list = ls())

library(egssimtools)
library(tidyverse)
library(ggsim)
library(patchwork)

# This is saving a hundred genomes

root <- "/media/raphael/bigass/simulations/EGS/genomes/"
data <- readRDS(paste0(root, "simulations.rds"))
backup <- data
data <- backup

head(data)

data <- data %>% filter(time == 19900)

ggplot
