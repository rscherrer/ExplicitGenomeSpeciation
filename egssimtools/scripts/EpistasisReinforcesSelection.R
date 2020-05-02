# We pick three random seeds that we simulate under ecological conditions such
# that competitive speciation occurs under additive and epistatic genetics

rm(list = ls())

library(egssimtools)
library(ggsim)
library(tidyverse)
library(cowplot)

# The pattern is the same with example1, 2 and 3

data <- list(
  additive = "/media/raphael/bigass/simulations/EGS/home/seed1/hsymmetry0/ecosel0/dispersal0.01/additive",
  epistatic = "/media/raphael/bigass/simulations/EGS/home/seed1/hsymmetry0/ecosel0/dispersal0.01/epistatic"
) %>%
  map_dfr(
    read_data,
    c("time", "individual_trait"),
    by = c(1, 3), dupl = list("population_size", 1), .id = "genetics")

head(data)

data %>%
  filter(time == 1990) %>%
  group_by(genetics) %>%
  mutate(individual_trait1 = individual_trait1 %>% scale(., center = TRUE, scale = FALSE)) %>%
  ungroup() %>%
  ggplot() +
  geom_density(aes(x = individual_trait1)) +
  facet_wrap(. ~ genetics)

data %>%
  ggplot() +
  geom_bin2d(aes(x = time, y = individual_trait1), bins = 50) +
  facet_wrap(. ~ genetics) +
  scale_fill_continuous(type = "viridis")

# Conclusion: epistasis seems to consistently constrain phenotypic variation in
# the nascent ecotypes.
