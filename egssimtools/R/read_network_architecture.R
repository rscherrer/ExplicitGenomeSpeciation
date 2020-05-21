#' Read gene network architecture
#'
#' Reads a genetic architecture and stores it into a `tbl_graph` object
#'
#' @param folder Path to the simulation
#' @param filename Optional architecture file name
#'
#' @return A `tbl_graph` object
#'
#' @export

read_network_architecture <- function(folder, filename = "architecture.txt") {

  library(tidyverse)
  library(tidygraph)

  # Read the nodes
  nodes <- read_genome_architecture(folder, filename)

  # Read the edges
  edges <- read_architecture(folder, filename)
  edges <- map_dfr(edges$networks, ~ do.call("tibble", .x), .id = "trait")
  edges <- edges %>% mutate(trait = as.numeric(factor(trait)) - 1)
  edges <- edges %>% rename(from = edges0, to = edges1, weight = weights)
  edges <- edges %>% mutate(from = from + 1, to = to + 1) # C++ indexing
  edges <- edges %>% select(from, to, weight, trait)

  # Return a tidygraph object suitable for plotting networks
  tbl_graph(nodes = nodes, edges = edges)

}
