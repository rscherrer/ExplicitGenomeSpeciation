# In this script we test the functions of the package

rm(list = ls())

library(egssimtools)
library(tidyverse)


data <- lapply(sprintf("/media/raphael/bigass/simulations/EGS/EGS_sim%s", seq(6)), function(root) {
  collect_simulations(
    root,
    c("EI", "SI", "RI"),
    parnames = c("ecosel", "hsymmetry", "dispersal", "mutation", "scaleA", "scaleI"),
    to_numeric = c("ecosel", "hsymmetry", "dispersal", "mutation")
  )
})
data <- do.call("rbind", data)

head(data)

# here is our heatmap figure, no need for a special function...
data %>%
  group_by(hsymmetry, ecosel, dispersal, mutation, scaleA, scaleI, simulation) %>%
  summarize(EI = last(EI), SI = last(SI), RI = last(RI)) %>%
  group_by(hsymmetry, ecosel, dispersal, mutation, scaleA, scaleI) %>%
  summarize(EI = mean(EI), SI = mean(SI), RI = mean(RI)) %>%
  ungroup() %>%
  ggplot(aes(x = hsymmetry, y = ecosel, fill = EI)) +
  facet_grid(mutation ~ scaleI) +
  geom_tile()

smr <- data %>%
  filter(scaleI == "1 0 0", mutation == 0.001) %>%
  mutate(ecosel = fct_rev(factor(ecosel)), hsymmetry = factor(hsymmetry)) %>%
  group_by(hsymmetry, ecosel, dispersal, mutation, scaleA, scaleI, simulation) %>%
  mutate(color = last(SI))

ggplot(smr, aes(x = time, y = EI, alpha = simulation, color = color)) %>%
  facettize(
    smr,
    facet_rows = "ecosel",
    facet_cols = "hsymmetry",
    label_facets = TRUE,
    facet_prefixes = c("s", "h")
  ) +
  geom_line() +
  scale_alpha_manual(values = runif(2000, min = 0.49, max = 0.51)) + # hack
  guides(alpha = FALSE) +
  scale_color_gradient(low = "black", high = "coral")




root <- sprintf("/media/raphael/bigass/simulations/EGS/EGS_sim%s", 2)
find_completed(root)

data <- collect_simulations(root, c("EI", "SI", "RI"), parnames = "ecosel")

status <- collect_status(root)
table(status)
find_extinct(root)

params <- collect_parameters(root, parnames = c("ecosel", "hsymmetry", "dispersal", "mutation", "scaleA", "scaleI"))
sims <- list.files(root, pattern = "sim_", full.names = TRUE)

status <- collect_status(root)
table(status)

get_status(list.files(root, pattern = "^sim_", full.names = TRUE)[1])

nrow(data)
nrow(params)

simulation <- "../build/release"

# Low-level functions: work on single simulations

read_paramfile(paste0(simulation, "/paramlog.txt"))
read_parameters(simulation)
read_binary(paste0(simulation, "/time.dat"))
read_data(simulation, "time")
read_time(simulation)
read_population_size(simulation)
read_population(simulation, "EI")
is_extinct(simulation)
is_missing(simulation)
read_archfile(paste0(simulation, "/architecture.txt"))
read_architecture(simulation)
read_individuals(simulation, "individual_trait")
read_ecotypes(simulation)
read_habitats(simulation)
read_loci(simulation, "genome_Fst")
read_edges(simulation, "network_corgen")
read_resources(simulation)
read_means(simulation)
read_ecotype_means(simulation)
read_ecotype_sizes(simulation)
plot_trait_density(simulation, trait = 1)
plot_trait_density2D(simulation, traits = c(2, 1), t = 100)


# Medium-level functions: can combine and process data from different simulations but the output needs to be processed to produce full-fledged figures and results

root <- "/media/raphael/bigass/simulations/EGS/EGS_sim1"

find_extinct(root, pattern = "sim_")
find_missing(root, pattern = "sim_")
simulations <- find_extant(root, pattern = "sim_", pb = FALSE)
collect_parameters(simulations[1:4], parnames = c("ecosel", "hsymmetry"))
collect_simulations(simulations[1:4], "RI", parnames = c("ecosel", "hsymmetry"))

data <- collect_simulations(simulations, c("RI", "EI", "SI"), parnames = c("ecosel", "hsymmetry"))
head(data)

# High-level functions: produce figures and overall results over multiple simulations


add_summaries <- function(data) data %>% group_by(simulation) %>% mutate(x = last(RI)) %>% ungroup() %>% select(x)


simulations <- find_extant(root)
collect_parameters(simulations)

data <- collect_simulations(simulations[1:4], c("RI", "EI", "SI"), parnames = c("ecosel", "hsymmetry"))
head(data)

params <- collect_parameters(simulations, c("ecosel", "hsymmetry"))
