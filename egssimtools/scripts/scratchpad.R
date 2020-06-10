rm(list = ls())

library(egssimtools)
library(tidyverse)
library(patchwork)
library(tidygraph)
library(ggraph)

root <- "data/example_1"

#### Example 1: plot ecological divergence through time ####

data <- read_sim(root, "EI")
ggplot(data, aes(x = time, y = EI)) +
  geom_line()

#### Example 2: pivoting to plot multiple speciation metrics ####

data <- read_sim(root, c("EI", "RI", "SI"))
data <- pivot_data(data, c("EI", "RI", "SI"))
ggplot(data, aes(x = time, y = value, color = variable)) +
  geom_line()

#### Example 3: splitting into several columns ####

data <- read_sim(root, "Fst", by = 3)
data <- pivot_data(data, paste0("Fst", 1:3), newnames = 0:2)
data <- data %>% rename(trait = "variable", Fst = "value")
ggplot(data, aes(x = time, y = Fst, color = trait)) +
  geom_line() +
  facet_grid(. ~ trait)

#### Example 4: individual-wise variables ####

data <- read_pop(root, "individual_trait", by = 3)
newnames <- paste0("trait ", 0:2)
data <- pivot_data(data, paste0("individual_trait", 1:3), newnames = newnames)
data <- data %>% rename(trait = "variable")
ggplot(data, aes(x = time, y = value)) +
  geom_bin2d(bins = 20) +
  facet_grid(. ~ trait)

#### Example 5: locus-wise variables ####

data <- read_genome(root, "genome_Fst")
data <- data %>% filter(time == last(time))
ggplot(data, aes(x = locus, y = genome_Fst)) +
  geom_point()

data <- read_genome(root, "genome_Fst")
ggplot(data, aes(x = time, y = locus, fill = genome_Fst)) +
  geom_tile()

#### Example 6: the genetic architecture ####

data <- read_genome(root, "genome_Fst", architecture = TRUE)
data <- data %>% filter(time == last(time))
ggplot(data, aes(x = location, y = genome_Fst, color = trait)) +
  geom_point()

#### Example 7: combining plots

data <- read_genome(root, "genome_Fst", architecture = TRUE)
data <- data %>% filter(time == last(time))
p1 <- ggplot(data, aes(x = location, y = genome_Fst, color = trait)) +
  geom_point()
p2 <- ggplot(data, aes(x = trait, y = genome_Fst, color = trait)) +
  geom_violin()
wrap_plots(p1, p2)

#### Example 8: more complex example with smoothing and alpha scales ####

data <- read_genome(root, "genome_Fst", architecture = TRUE)
data <- smoothen_data(data, "time", "genome_Fst", span = 0.3, line = "locus")
ggplot(data, aes(x = time, y = genome_Fst, alpha = factor(locus), color = trait)) +
  geom_line() +
  facet_grid(. ~ trait) +
  guides(alpha = FALSE)

### Example 9: edge-wise data ####

data <- read_network(root, "network_corfreq", architecture = TRUE)
ggplot(data, aes(x = trait, y = network_corfreq, color = trait)) +
  geom_violin()

#### Example 10: the reading functions, what happens in the background ####

read_data(
  root,
  c("time", "genome_Fst", "genome_Cst"),
  by = c(1, 3, 1),
  dupl = c(300, 1, 1)
)

read_genome(root, c("genome_Fst", "genome_Cst"))

read_data(
  root,
  c("time", "individual_trait", "individual_ecotype"),
  by = c(1, 3, 1),
  dupl = list("population_size", 1, 1)
)

read_pop(root, paste0("individual_", c("trait", "ecotype")), by = c(3, 1))

#### Example 11: plot a gene network ####

arch <- read_arch_network(root)
data <- read_genome(root, "genome_Fst", architecture = TRUE)
data <- data %>% filter(time == last(time))
arch <- arch %>% activate(nodes) %>% right_join(data)
ggraph(arch, layout = "graphopt", charge = 0.1, mass = 30, niter = 100000) +
  geom_edge_link(aes(color = trait), width = 2, alpha = 0.6) +
  geom_node_point(aes(fill = trait, size = genome_Fst), shape = 21) +
  scale_size_continuous(range = c(2, 7)) +
  scale_alpha(range = c(0.6, 1)) +
  labs(size = parse(text = "F[ST]"), fill = "trait") +
  theme_void() +
  guides(edge_color = FALSE)

#### Example 12: combining simulations ####

root <- "data"
data <- combine_data(
  root, pattern = "example", level = 1, type = "genome",
  variables = "genome_Fst", architecture = TRUE,
  parnames = c("ecosel", "seed")
)
data <- data %>% filter(time == last(time))
data <- data %>% mutate(sim = str_replace(sim, "^", "simulation "))
ggplot(data, aes(x = trait, y = genome_Fst)) +
  geom_violin() +
  facet_grid(sim ~ ecosel)

#### Example 13: the split-apply-combine routine ####

root <- "data"
data <- combine_data(
  root, pattern = "example", level = 1, type = "genome",
  variables = "genome_Fst", architecture = TRUE
)

plot_this <- function(data) {

  ggplot(
    data,
    aes(x = time, y = genome_Fst, color = trait, alpha = factor(locus))
  ) +
    geom_line() +
    facet_grid(. ~ trait) +
    guides(alpha = FALSE)

}

data <- data %>%
  group_by(sim) %>%
  nest() %>%
  mutate(fig = map(data, plot_this))
data$fig[[1]]

# Save each figure separately
data <- data %>% mutate(figname = sprintf("sim%s.png", sim))
save_this <- function(x, y) ggsave(x, y, width = 4, height = 3, dpi = 300)
data %>% mutate(saved = walk2(figname, fig, save_this))

# Or combine into one
wrap_plots(data$fig) +
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'A')
