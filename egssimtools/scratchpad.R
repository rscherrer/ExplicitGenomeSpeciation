rm(list = ls())

library(egssimtools)
library(tidyverse)
library(ggsim)
library(patchwork)

root <- "data"

# Fetch simulations folders
fetch_dirs(roots = root, level = 1)



data <- readRDS(paste0(root, "simulations.rds"))
backup <- data
data <- backup

head(data)

data <- data %>% filter(time == 19900)

ggplot
