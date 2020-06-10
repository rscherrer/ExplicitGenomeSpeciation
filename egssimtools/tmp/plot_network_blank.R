#' Plot gene network without evolving variable
#'
#' @param root Path to the simulation folder
#' @param y Optional name of the variable to map as node size aesthetics.
#' Should be in the fields returned by `read_genome_architecture`.
#' @param archfile Optional architecture file name
#'
#' @return A ggplot
#'
#' @examples
#'
#' root <- system.file("extdata", "example_1", package = "egssimtools")
#' plot_network_blank(root, y = "degree")
#'
#' @export

plot_network_blank <- function(root, y = NULL, archfile = "architecture.txt") {

  arch <- read_network_architecture(root)
  arch <- arch %>%
    tidygraph::activate(edges) %>%
    dplyr::mutate(trait = factor(trait))

  p <- ggraph::ggraph(
    arch, layout = "graphopt", charge = 0.1, mass = 30, niter = 100000
  ) +
    ggraph::geom_edge_link(
      ggplot2::aes(color = factor(trait)), width = 2, alpha = 0.6
    ) +
    ggraph::geom_node_point(
      ggplot2::aes_string(fill = "trait", size = y), shape = 21
    ) +
    ggplot2::scale_size_continuous(range = c(2, 7)) +
    ggplot2::scale_alpha(range = c(0.6, 1)) +
    ggplot2::labs(size = y, fill = "trait") +
    ggplot2::theme_void() +
    ggplot2::guides(edge_color = FALSE)

  p

}
