rm(list = ls())

library(egssimtools)
library(tidyverse)
library(patchwork)

roots <- fetch_dirs("data", pattern = "example", level = 1)
root <- roots[1]

#### Example 1: plot speciation metrics through time ####

data <- read_sim(root, c("EI", "RI", "SI"))
data <- pivot_data(data, c("EI", "RI", "SI"))

ggplot(data, aes(x = time, y = value, color = variable)) +
  geom_line()

#### Example 2: plot trait distributions through time ####

data <- read_pop(root, "individual_trait", by = 3)
cols <- paste0("individual_trait", 1:3)
newnames <- 0:2
data <- pivot_data(data, cols, newnames)
data <- data %>% rename(trait = "variable")

ggplot(data, aes(x = time, y = value)) +
  geom_bin2d(bins = 20) +
  facet_grid(trait ~ .)

#### Example 3: compare simulation metrics across simulations ####

variables <- c("time", "EI", "RI", "SI")
data <- collect_sims(roots, variables, check_extant = FALSE, level = 0)
data <- pivot_data(data, variables[-1])

ggplot(data, aes(x = time, y = value, color = sim)) +
  geom_line() +
  facet_grid(. ~ variable)

#### Example 4: compare genome scans across simulations ####

data <- collect_sims(
  roots, c("time", "genome_Fst"), dupl = c(300, 1), check_extant = FALSE,
  level = 0, architecture = TRUE
)
data <- data %>% filter(time == last(time))

ggplot(data, aes(x = locus, y = genome_Fst, color = trait)) +
  geom_point() +
  facet_grid(sim ~ .)

#### Example 5: compare Fst through time across traits and simulations ####

data <- collect_sims(
  roots, c("time", "genome_Fst"), dupl = c(300, 1), check_extant = FALSE,
  level = 0, architecture = TRUE
)

plot_this <- function(data) {

  ggplot(data, aes(x = time, y = genome_Fst, alpha = factor(locus))) +
    geom_line() +
    guides(alpha = FALSE) +
    facet_grid(trait ~ .)

}

data <- data %>%
  group_by(sim) %>%
  nest() %>%
  mutate(fig = map(data, plot_this))

wrap_plots(data$fig)

#### Example 6: Fst heatmap through time ####

data <- read_genome(root, "genome_Fst")
ggplot(data, aes(x = time, y = locus, fill = genome_Fst)) +
  geom_tile()

