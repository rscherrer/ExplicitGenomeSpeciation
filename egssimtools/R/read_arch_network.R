#' Read gene network architecture
#'
#' Reads a genetic architecture and stores it into a `tbl_graph` object
#'
#' @param folder Path to the simulation
#' @param filename Optional architecture file name
#' @param as_list Whether to return the output at a list of two tibbles
#' instead of a `tbl_graph`.
#'
#' @return A `tbl_graph` object (useful for plotting using the `ggraph`
#' package), or a list of data frames if applicable
#'
#' @examples
#'
#' \dontrun{
#'
#' # Location of the simulation folder
#' root <- "data/example_1"
#'
#' read_arch_network(root)
#'
#' }
#'
#' @export

read_arch_network <- function(
  folder,
  filename = "architecture.txt",
  as_list = FALSE
) {

  # Read the nodes
  nodes <- read_arch_genome(folder, filename)

  # Read the edges
  edges <- read_arch(folder, filename)
  edges <- purrr::map_dfr(
    edges$networks, ~ do.call("data.frame", .x), .id = "trait"
  )
  edges <- edges %>% dplyr::mutate(trait = factor(as.numeric(factor(trait))-1))
  edges <- edges %>% dplyr::rename(from = edges0, to = edges1, weight = weights)
  edges <- edges %>% dplyr::mutate(from = from + 1, to = to + 1) # C++ indexing
  edges <- edges %>% dplyr::mutate(edge = seq(nrow(edges)))
  edges <- edges %>% dplyr::select(from, to, weight, trait, edge)

  # Return a tidygraph object suitable for plotting networks
  network <- tidygraph::tbl_graph(nodes = nodes, edges = edges)
  if (as_list) network <- as.list(network)
  return (network)

}
