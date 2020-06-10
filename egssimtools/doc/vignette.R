## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message = FALSE---------------------------------------------------------
library(egssimtools)
library(tidyverse)
library(patchwork)

## -----------------------------------------------------------------------------
roots <- fetch_dirs("../data", pattern = "example", level = 1)
roots
# we are within the "vignettes" folder, hence the ".."

## -----------------------------------------------------------------------------
root <- roots[1]

## -----------------------------------------------------------------------------
data <- read_sim(root, c("EI", "RI", "SI"))
data

## -----------------------------------------------------------------------------
data <- pivot_data(data, c("EI", "RI", "SI"))
data

## ---- fig.width = 3, fig.height = 2, fig.align = "center"---------------------
ggplot(data, aes(x = time, y = value, color = variable)) +
  geom_line()

## -----------------------------------------------------------------------------
data <- read_pop(root, "individual_trait", by = 3)
data

## -----------------------------------------------------------------------------
cols <- paste0("individual_trait", 1:3)
newnames <- paste0("trait ", 0:2) # to match the C++ numbering of traits
data <- pivot_data(data, cols, newnames = newnames)
data <- data %>% rename(trait = "variable")
data

## ---- fig.width = 5, fig.height = 2, fig.align = "center"---------------------
ggplot(data, aes(x = time, y = value)) +
  geom_bin2d(bins = 20) +
  facet_grid(. ~ trait)

## -----------------------------------------------------------------------------
variables <- c("time", "Fst")
data <- collect_data(
  roots, variables, by = c(1, 3), check_extant = FALSE, level = 0
)
data

## ---- fig.width = 5, fig.height = 2, fig.align = "center"---------------------
data <- pivot_data(data, paste0("Fst", 1:3), newnames = newnames)
data <- data %>% rename(trait = "variable")
data
ggplot(data, aes(x = time, y = value, color = sim)) +
  geom_line() +
  facet_grid(. ~ trait) +
  ylab(parse(text = "F[ST]"))

## -----------------------------------------------------------------------------
data <- collect_data(
  roots, c("time", "genome_Fst"), dupl = c(300, 1), check_extant = FALSE,
  level = 0, architecture = TRUE
)
data

## ---- fig.width = 5, fig.height = 3, fig.align = "center"---------------------
data <- data %>% filter(time == last(time))
ggplot(data, aes(x = locus, y = genome_Fst, color = trait)) +
  geom_point() +
  facet_grid(sim ~ .)

## ---- fig.width = 3, fig.height = 3-------------------------------------------
data <- read_genome(root, "genome_Fst")
ggplot(data, aes(x = time, y = locus, fill = genome_Fst)) +
  geom_tile()

## -----------------------------------------------------------------------------
data <- collect_data(
  roots, c("time", "genome_Fst"), dupl = c(300, 1), check_extant = FALSE,
  level = 0, architecture = TRUE
)
data <- data %>% mutate(trait = str_replace(trait, "^", "trait "))
data

## -----------------------------------------------------------------------------
plot_this <- function(data) {

  ggplot(data, aes(x = time, y = genome_Fst, alpha = factor(locus))) +
    geom_line() +
    guides(alpha = FALSE) +
    facet_grid(trait ~ .)

}

## -----------------------------------------------------------------------------
data <- data %>%
  group_by(sim) %>%
  nest() %>%
  mutate(fig = map(data, plot_this))
data

## ---- warning = FALSE, fig.width = 3, fig.height = 3, fig.align = "center"----
data$fig[[1]]

## ---- warning = FALSE, fig.width = 7, fig.height = 3, fig.align = "center"----
wrap_plots(data$fig)

## -----------------------------------------------------------------------------
data$figname <- sprintf("sim%s.png", 1:3)
save_this <- function(figname, fig) {
  ggsave(figname, fig, width = 4, height = 3, dpi = 300)
}
#data <- data %>% mutate(saved = walk2(figname, fig, save_this))
# uncomment to actually save the plots

## -----------------------------------------------------------------------------
data <- read_data(
  root, 
  variables = c("time", "individual_trait", "individual_ecotype"), 
  by = c(1, 3, 1),
  dupl = list("population_size", 1, 1),
  parnames = c("ecosel", "hsymmetry")
)
data

## -----------------------------------------------------------------------------
data <- read_pop(
  root, 
  variables = c("individual_trait", "individual_ecotype"),
  by = c(3, 1),
  parnames = c("ecosel", "hsymmetry")
)
data

