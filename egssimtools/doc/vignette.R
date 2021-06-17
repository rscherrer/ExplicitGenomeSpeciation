## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#"
)

## ---- message = FALSE---------------------------------------------------------
library(egssimtools)
library(tidyverse)
library(patchwork)
library(tidygraph)
library(ggraph)

## -----------------------------------------------------------------------------
root <- "../data/example_1"

## -----------------------------------------------------------------------------
data <- read_sim(root, "EI")
data

## ---- fig.width = 3, fig.height = 2, fig.align = "center"---------------------
ggplot(data, aes(x = time, y = EI)) +
  geom_line()

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
data <- read_sim(root, "Fst", by = 3)
data

## -----------------------------------------------------------------------------
data <- pivot_data(data, paste0("Fst", 1:3), newnames = paste0("trait ", 0:2))
data <- data %>% rename(trait = "variable", Fst = "value")
data

## ---- fig.width = 6, fig.height = 2, fig.align = "center"---------------------
ggplot(data, aes(x = time, y = Fst, color = trait)) +
  geom_line() +
  facet_grid(. ~ trait) +
  theme(legend.position = "none")

## -----------------------------------------------------------------------------
data <- read_pop(root, "individual_trait", by = 3)
data

## ---- fig.width = 6, fig.height = 2, fig.align = "center"---------------------
newnames <- paste0("trait ", 0:2)
data <- pivot_data(data, paste0("individual_trait", 1:3), newnames = newnames)
data <- data %>% rename(trait = "variable")
ggplot(data, aes(x = time, y = value)) +
  geom_bin2d(bins = 20) +
  facet_grid(. ~ trait)

## -----------------------------------------------------------------------------
data <- read_genome(root, "genome_Fst")
data

## -----------------------------------------------------------------------------
read_sim(root, "genome_Fst", by = 300)

## ---- fig.width = 4, fig.height = 2, fig.align = "center"---------------------
data <- data %>% filter(time == last(time))
ggplot(data, aes(x = locus, y = genome_Fst)) +
  geom_point()

## ---- fig.width = 4, fig.height = 4, fig.align = "center"---------------------
data <- read_genome(root, "genome_Fst")
ggplot(data, aes(x = time, y = locus, fill = genome_Fst)) +
  geom_tile()

## -----------------------------------------------------------------------------
data <- read_genome(root, "genome_Fst", architecture = TRUE)
data

## ---- fig.width = 4, fig.height = 2, fig.align = "center"---------------------
data <- data %>% filter(time == last(time))
ggplot(data, aes(x = location, y = genome_Fst, color = trait)) +
  geom_point()

## ---- fig.width = 5, fig.height = 2, fig.align = "center"---------------------
data <- read_genome(root, "genome_Fst", architecture = TRUE)
data <- data %>% filter(time == last(time))
p1 <- ggplot(data, aes(x = location, y = genome_Fst, color = trait)) +
  geom_point() +
  theme(legend.position = "none")
p2 <- ggplot(data, aes(x = trait, y = genome_Fst, color = trait)) +
  geom_violin()

# From patchwork
wrap_plots(p1, p2, widths = c(2, 1)) + plot_annotation(tag_levels = "A")

## ---- fig.width = 6, fig.height = 2, fig.align = "center"---------------------
data <- read_genome(root, "genome_Fst", architecture = TRUE)
data <- smoothen_data(data, "time", "genome_Fst", span = 0.3, line = "locus")
ggplot(data, aes(x = time, y = genome_Fst, alpha = factor(locus), color = trait)) +
  geom_line() +
  facet_grid(. ~ trait) +
  scale_alpha_manual(values = runif(length(unique(data$locus)), 0.79, 0.81)) +
  guides(alpha = FALSE)

## -----------------------------------------------------------------------------
data <- read_network(root, "network_corfreq", architecture = TRUE)
data

## ---- fig.width = 3, fig.height = 2, fig.align = "center"---------------------
ggplot(data, aes(x = trait, y = network_corfreq, color = trait)) +
  geom_violin()

## -----------------------------------------------------------------------------
read_data(
  root,
  c("time", "genome_Fst", "genome_Cst"),
  by = c(1, 1, 1),
  dupl = c(300, 1, 1)
)

## -----------------------------------------------------------------------------
read_genome(root, c("genome_Fst", "genome_Cst"))

## -----------------------------------------------------------------------------
read_data(
  root,
  c("time", "individual_trait", "individual_ecotype"),
  by = c(1, 3, 1),
  dupl = list("population_size", 1, 1)
)

## -----------------------------------------------------------------------------
read_pop(root, paste0("individual_", c("trait", "ecotype")), by = c(3, 1))

## -----------------------------------------------------------------------------
arch <- read_arch_network(root)
arch

## ---- message = FALSE---------------------------------------------------------
data_n <- read_genome(root, "genome_Fst", architecture = TRUE)
data_n <- data_n %>% filter(time == last(time))

## ---- message = FALSE---------------------------------------------------------
arch <- arch %>% activate(nodes) %>% right_join(data_n)
arch

## ---- message = FALSE---------------------------------------------------------
data_e <- read_network(root, "network_corfreq", architecture = TRUE)
data_e <- data_e %>% filter(time == last(time))
arch <- arch %>% activate(edges) %>% right_join(data_e)

## ----  fig.width = 4, fig.height = 4, fig.align = "center"--------------------
# This may take a while
ggraph(arch, layout = "graphopt", charge = 0.1, mass = 30, niter = 100000) +
  geom_edge_link(aes(color = trait), width = 2, alpha = 0.6) +
  geom_node_point(aes(fill = trait, size = genome_Fst), shape = 21) +
  scale_size_continuous(range = c(2, 7)) +
  scale_alpha(range = c(0.6, 1)) +
  labs(size = parse(text = "F[ST]"), fill = "trait") +
  theme_void() +
  guides(edge_color = FALSE)

## ---- message = FALSE---------------------------------------------------------
# First we reset the working directory to where multiple simulations are
root <- "../data"

data <- combine_data(
  root, pattern = "example", level = 1, type = "genome",
  variables = "genome_Fst", architecture = TRUE,
  parnames = c("ecosel", "hsymmetry")
)
data

## ---- fig.width = 4, fig.height = 4, fig.align = "center"---------------------
data <- data %>% filter(time == last(time))
data <- data %>% mutate(sim = str_replace(sim, "^", "simulation "))
data <- data %>% mutate(ecosel = str_replace(ecosel, "^", "ecosel = "))
ggplot(data, aes(x = trait, y = genome_Fst)) +
  geom_violin() +
  facet_grid(sim ~ ecosel)

## ---- message = FALSE---------------------------------------------------------
data <- combine_data(
  root, pattern = "example", level = 1, type = "genome",
  variables = "genome_Fst", architecture = TRUE
)

## -----------------------------------------------------------------------------
data <- data %>%
  group_by(sim) %>%
  nest()
data

## -----------------------------------------------------------------------------
plot_this <- function(data) {

  ggplot(
    data,
    aes(x = time, y = genome_Fst, color = trait, alpha = factor(locus))
  ); geom_line(); facet_grid(. ~ trait); guides(alpha = FALSE)

}

## -----------------------------------------------------------------------------
data <- data %>% mutate(fig = map(data, plot_this))
data

## ---- warning = FALSE, fig.width = 6, fig.height = 2, fig.align = "center"----
data$fig[[1]]

## -----------------------------------------------------------------------------
data <- data %>% mutate(figname = sprintf("sim%s.png", sim))
data
save_this <- function(x, y) ggsave(x, y, width = 4, height = 3, dpi = 300)

# This is commented to not save when rendering the vignette
# data %>% mutate(saved = walk2(figname, fig, save_this))

## -----------------------------------------------------------------------------
roots <- fetch_dirs("../data", pattern = "example_", level = 1)
roots

## -----------------------------------------------------------------------------
plot_this <- function(root) {
  
  # This is the only line added
  data <- read_genome(root, "genome_Fst", architecture = TRUE)

  ggplot(
    data,
    aes(x = time, y = genome_Fst, color = trait, alpha = factor(locus))
  ); geom_line(); facet_grid(. ~ trait); guides(alpha = FALSE)

}

## ---- message = FALSE---------------------------------------------------------
plots <- map(roots, plot_this)

## ---- warning = FALSE, fig.width = 6, fig.height = 2, fig.align = "center"----
plots[[1]]

## ---- warning = FALSE, fig.width = 6, fig.height = 5, fig.align = "center"----
# We remove the legend of all plots but one
for (i in 2:length(plots)) {
  plots[[i]] <- plots[[i]] + theme(legend.position = "none")
}

# Then we assemble the plots with a common legend
wrap_plots(plots) +
  plot_layout(guides = 'collect', nrow = length(plots)) +
  plot_annotation(tag_levels = 'A')

