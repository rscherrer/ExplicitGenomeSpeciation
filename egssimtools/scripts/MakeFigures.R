# Characterization of the model behavior throughout parameter space
# Use this script to generate simple, combined or facetted figures
# showing the model behavior across parameter space.

# The script is organized in three parts:
# 1. Collect simulation data and make a master database
# 2. Plot heatmaps throughout parameter space
# 3. Plot simulation trajectories for specific subparts of parameter space

rm(list = ls())

library(egssimtools)
library(tidyverse)
library(cowplot)
library(ggsim)

# First choose your settings

make_database <- FALSE
database <- "data/simulations.rds"
variables <- c("EI", "RI", "SI")
colors <- c("lightgreen", "lightblue", "coral")
ylabs <- c("Ecological divergence", "Reproductive isolation", "Spatial isolation")

#### 1. Make database ####

# This chunk will loop through the simulation folders defined in root
# and extract the variables and parameters of interest into a data frame
# If you make_database is FALSE, it will read a previously made database from
# the address you specified

"/media/raphael/bigass/simulations/EGS/EGS_sim1" %>%
  fetch_dirs("sim_", level = 1) %>%
  .[1:5] %>%
  map(read_data, c("EI"))

if (make_database) {

  root <- "/media/raphael/bigass/simulations/EGS"
  data <- collect_simulations(
    root, level = 2, variables = variables,
    parnames = c("ecosel", "hsymmetry", "dispersal", "mutation", "scaleA", "scaleI"),
    to_numeric = c("ecosel", "hsymmetry", "dispersal", "mutation"),
    as_address = TRUE
  )
  saveRDS(data, database)

} else data <- readRDS(database)

#### 2. Plot heatmaps ####

# Use this chunk to make heatmaps of the simulation outcomes across parameter
# space.

list(variables, colors) %>%
  pmap(
    function(variable, color) {
      data %>%
        mutate(
          scaleI = scaleI %>% str_replace(" .*$", "") %>% str_replace("^", "sigma[I]=="),
          mutation = mutation %>% factor %>% str_replace("^", 'mu=="') %>% str_replace("$", '"') %>% str_replace("1e-04", "0.0001") %>% fct_rev
        ) %>%
        ggheatmap(variable, "hsymmetry", "ecosel", "simulation", c(last, function(x) length(which(x > 0.9))), keep = c("mutation", "scaleI")) +
        facet_grid(mutation ~ scaleI, labeller = label_parsed) +
        scale_fill_gradient(low = "black", high = color) +
        theme_bw() +
        labs(x = "Habitat symmetry", y = "Divergent selection", fill = variable)
    }
  ) %>%
  plot_grid(plotlist = ., nrow = 3, labels = c("A", "B", "C"))
#ggsave("figures/heatmaps.png", width = 4, height = 7, dpi = 300)

#### 3. Plot trajectories ####

# Use this chunk to plot trajectories of the simulations

data %>%
  mutate(
    ecosel = ecosel %>% factor %>% str_replace("^", "s = ") %>% fct_rev,
    hsymmetry = hsymmetry %>% factor %>% str_replace("^", "h = "),
    time = time / 1000
  ) %>%
  group_by(mutation, scaleI) %>%
  group_map(function(df, ...) {

    filename <- "figures/trajectories_mutation_%s_scaleI_%s.png" %>%
      sprintf(df$mutation[1], df$scaleI[1])

    print(filename)

    list(variables, colors, ylabs) %>%
      pmap(function(variable, color, label) {

        df %>%
          group_by(simulation) %>%
          mutate(last = last(get(variable))) %>%
          ungroup() %>%
          gglineplot(x = "time", y = variable, line = "simulation") +
          aes(color = last) +
          scale_color_gradient(low = "black", high = color) +
          guides(color = FALSE) +
          labs(x = "Time (\U02715 1000 generations)", y = label) +
          facet_grid(ecosel ~ hsymmetry)

      }) %>%
      plot_grid(plotlist = ., ncol = 3, labels = c("A", "B", "C")) %>%
      ggsave(filename, ., width = 13, height = 9, dpi = 400)
  }, keep = TRUE)
