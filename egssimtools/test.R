# In this script we test the functions of the package

rm(list = ls())

library(egssimtools)
library(tidyverse)
library(cowplot)

simulation <- "../build/release"

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
find_extinct("../build", pattern = "release")
find_missing("../build", pattern = "release")
read_resources(simulation)
read_means(simulation)
read_ecotype_means(simulation)
read_ecotype_sizes(simulation)

# Plot a simulation through time
data <- read_individuals(simulation, "individual_trait")
data <- lapply(data, function(data) sapply(data, "[", 1))
data <- data.frame(time = mrep(as.numeric(names(data)), sapply(data, length)), data = do.call("c", data))

p <- ggplot(data, aes(x = time, y = data)) +
  geom_bin2d(bins = 15) +
  theme_bw() +
  scale_fill_continuous(type = "viridis") +
  xlab("Time (generations)") +
  ylab("Ecological trait") +
  labs(fill = "Count")
p

matepref <- read_individuals(simulation, "individual_trait")
matepref <- do.call("c", lapply(matepref, function(data) sapply(data, "[", 2)))
data$matepref <- matepref

p2 <- ggplot(data %>% filter(time == 2900), aes(x = matepref, y = data)) +
  geom_bin2d(bins = 15) +
  theme_bw() +
  scale_fill_continuous(type = "viridis") +
  xlab("Mate preference trait") +
  ylab("Ecological trait") +
  labs(fill = "Count") +
  xlim(c(-3, 3))
p2

ylim <- c(-3, 3)
p <- p + theme(legend.position = "none") + ylim(ylim)
p2 <- p2 + ylab(NULL) + ylim(ylim)

p2 <- plot_grid(ggplot() + geom_blank() + theme_void(), p2, rel_widths = c(1, 11))

p3 <- plot_grid(p, p2, labels = c("A", "B"), rel_widths = c(1, 1))
p3

ggsave("../pics/example_speciation.png", p3, height = 2, width = 6, dpi = 300)
