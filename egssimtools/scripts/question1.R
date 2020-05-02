rm(list = ls())

library(egssimtools)
library(tidyverse)

# We see that when the mutation rate is high, reproductive isolation in asymmetric
# habitats evolves even less than when mutation is high. Why is that? To answer
# we compare the trait distributions through time in high and low mutation scenarios...
# Hypothesis: the high mutation rate increases the variance in ecological trait
# such that ecotypes are not ecologically isolated enough for RI to be advantageous.
# If this is true, we should see RI evolving at higher mutation rates, where ecotypes
# become more isolated.

root <- "/media/raphael/bigass/simulations/EGS_extra/EGS_sim3"

data <- collect_simulations(
  root, variables = c("EI", "RI", "SI"), parnames = c("ecosel", "hsymmetry"),
  to_numeric = c("ecosel", "hsymmetry")
)

head(data)

data %>%
  group_by(hsymmetry, ecosel, simulation) %>%
  summarize(EI = last(EI), RI = last(RI), SI = last(SI)) %>%
  ggplot(aes(x = hsymmetry, y = ecosel, fill = RI)) +
  geom_tile()

# EI still happens all over the place at high selection, still without RI, probably because
# selection is not strong enough. There is a 10-fold difference between the mutation rates
# so maye we need a 10-fold difference in selection intensity...? It would be nice to explore
# this with individual simulations before commiting to simulating 1000 more simulations
# on the cluster...

data <- readRDS("data/population_wide_data.rds")

head(data)

data %>%
  filter(hsymmetry == 0, scaleI == "0 0 0", ecosel == 1.2) %>%
  group_by(simulation, mutation) %>%
  summarize(EI = last(EI), RI = last(RI)) %>%
  split(f = .$mutation) %>% map(function(x) x[1, ]) %>%
  bind_rows()

# Let us rerun those simulations with the exact same seed but this time saving individual trait values
# and population sizes too...
# Rerun with the executable compiled  on the cluster (to be sure the results are exactly the same)
# Done!

sims <- smr$simulation

read_individuals("../build/release", 1)

plot_trait_density(sims[1], 2)
plot_trait_density(sims[2], 2)

plot_trait_histogram(sims[1], 1, t = 1000)
plot_trait_histogram(sims[2], 1, t = 1000)
