# Read the simulations and assemble them into a data set

rm(list = ls())

library(egssimtools)
library(tidyverse)
library(ggsim)
library(cowplot)
library(patchwork)

root <- "/media/raphael/bigass/simulations/EGS/genomes/"

# Read genome-wide variables

variables <- c("time", "Fst", "Gst", "Qst", "Cst", "varG", "varA", "varN", "varT")
by <- c(1, 3, 3, 3, 3, 3, 3, 3, 3)
dupl <- c(1, 1, 1, 1, 1, 1, 1, 1, 1)
parnames <- c("hsymmetry", "ecosel", "scaleI", "seed", "archfile")

data <- collect_sims(
  root, variables, by = by, dupl = dupl, parnames = parnames,
  check_extant = FALSE, level = 2, pattern = "sim"
)

# Check the data
head(data)

# Plot each genome-wide variable through time
# Facets correspond to different genetic / ecological scenarios
# Lines are different architectures and seeds
# Multiple plots are assembled per figure, one plot per trait
# One figure per variable



# Plot the data through time across ecological and genetic scenarios
plot_this <- function(y, title) {

  variable <- gsub("[0-9]", "", y)

  data %>%
    mutate(time = time / 1000) %>%
    gglineplot(x = "time", y = y, line = "sim") %>%
    facettize(
      "hsymmetry",
      "scaleI1",
      prepend = c("h = ", "sigma[I]=="),
      parsed = "scaleI1"
    ) +
    aes(color = archfile) +
    scale_color_brewer(palette = "Dark2") +
    theme(legend.position = "none") +
    xlab("Time (\U02715 1,000 generations)") +
    ylab(paste("Genome-wide", variable)) +
    ggtitle(title)

}

# Apply the plotting function to each trait to make a figure
figure_this <- function(variable) {

  y <- paste0(variable, 1:3)
  title <- c("A) Ecological trait", "B) Mating preference", "C) Neutral trait")
  p <- y %>% map2(title, plot_this)
  p[[1]] / p[[2]] / p[[3]]

}

# Now apply the figure-making function to all variables
figs <- variables[-1] %>% map(figure_this)

# Save the figures
fignames <- sprintf("figure_overview_%s.png", variables[-1])
figs %>% map2(fignames, ~ ggsave(.y, .x, width = 5, height = 10, dpi = 400))
