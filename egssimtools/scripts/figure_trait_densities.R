# Plot trait densities throughout the simulation

rm(list = ls())

library(egssimtools)
library(tidyverse)
library(cowplot)

simulation <- "../build/release"

# Ecological trait through time
p <- plot_trait_density(simulation, trait = 1) +
  ylab("Ecological trait")

# Ecological versus mating preference traits at the last generation
p2 <- plot_trait_density2D(simulation, traits = c(2, 1), t = 2900) +
  xlab("Mating preference") +
  ylab("Ecological trait") +
  xlim(c(-3, 3))

# Rearrange into a single figure
p <- p + theme(legend.position = "none") + ylim(c(-3, 3))
p2 <- p2 + ylab(NULL) + ylim(c(-3, 3))

# Blank plot to fill up space
p2. <- plot_grid(ggplot() + geom_blank() + theme_void(), p2, rel_widths = c(1, 11))

# Assemble the plots
p3 <- plot_grid(p, p2., labels = c("A", "B"), rel_widths = c(1, 1))
p3

ggsave("../pics/example_speciation.png", p3, height = 2, width = 6, dpi = 300)

####

root <- "/media/raphael/bigass/simulations/EGS/EGS_sim1"
simulations <- find_extant(root)
pars <- collect_parameters(simulations, c("ecosel", "hsymmetry"), to_numeric = c("ecosel", "hsymmetry"))
simulations <- simulations[pars$ecosel == 1.6 & pars$hsymmetry == 0]
plot_trait_density(simulations[1], "EI")

# Can I launch a simulation from R? but the seed on my machine won't be the same as the one on the cluster...
system("../build/release/EGS")

