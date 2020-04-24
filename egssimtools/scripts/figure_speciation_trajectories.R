rm(list = ls())

library(egssimtools)
library(tidyverse)
library(cowplot)

# Settings
variables <- c("EI", "RI", "SI")
colors <- c("lightgreen", "lightblue", "coral")

# Read the data
data <- readRDS("data/population_wide_data.rds")

# Parameter values
mutations <- unique(data$mutation)
scaleIs <- levels(data$scaleI)

# Variables to plot
variables <- c("EI", "RI", "SI")
colors <- c("lightgreen", "lightblue", "coral")
ylabs <- c("Ecological divergence", "Reproductive isolation", "Spatial isolation")

# For each combination of parameters...
lapply(mutations, function(mutation_value) {
  lapply(scaleIs, function(scaleI_value) {

    p <- mapply(function(variable, col, ylab) {

      # Summarize the data
      smr <- data %>%
        mutate(time = time / 1000) %>%
        filter(scaleI == scaleI_value, mutation == mutation_value) %>%
        mutate(ecosel = fct_rev(factor(ecosel)), hsymmetry = factor(hsymmetry)) %>%
        group_by(hsymmetry, ecosel, dispersal, mutation, scaleA, scaleI, simulation) %>%
        mutate(color = last(get(variable)))

      # Plot the trajectories
      p <- ggplot(smr, aes(x = time, y = get(variable), alpha = simulation, color = color)) %>%
        facettize(
          smr,
          facet_rows = "ecosel",
          facet_cols = "hsymmetry",
          label_facets = TRUE,
          facet_prefixes = c("s", "h")
        ) +
        geom_line() +
        xlab("Time (\U02715 1000 generations)") +
        ylab(ylab) +
        scale_x_continuous(n.breaks = 4) +
        scale_y_continuous(n.breaks = 4) +
        scale_alpha_manual(values = runif(2000, min = 0.49, max = 0.51)) + # hack
        guides(alpha = FALSE, color = FALSE) +
        scale_color_gradient(low = "black", high = col)

      return (p)

    }, variables, colors, ylabs, SIMPLIFY = FALSE)

    #p[[1]]
    #p[[2]]
    #p[[3]]

    # Combine and save the plots
    fig <- plot_grid(plotlist = p, ncol = 3, labels = c("A", "B", "C"))
    #fig
    ggsave(sprintf("figures/trajectories_mutation_%s_scaleI_%s.png", mutation_value, gsub(" ", "_", scaleI_value)), width = 12, height = 6, dpi = 300)

  })
})

