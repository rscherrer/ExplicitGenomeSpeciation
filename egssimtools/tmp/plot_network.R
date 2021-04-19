#' Plot gene network
#'
#' @param root Path to the simulation folder
#' @param y Name of the variable to map as node size aesthetics
#' @param t What time point? Defaults to last time point
#' @param archfile Optional architecture file name
#'
#' @return A ggplot
#'
#' @examples
#'
#' root <- system.file("extdata", "example_1", package = "egssimtools")
#' plot_network(root, "genome_Fst")
#'
#' @export

plot_network <- function(root, y, t = NULL, archfile = "architecture.txt") {

  data <- read_loci(root, y, architecture = TRUE, archfile = archfile)

  if (is.null(t)) t <- dplyr::last(data$time)
  data <- data %>% dplyr::filter(time == t)

  arch <- read_network_architecture(root)
  arch <- arch %>%
    tidygraph::activate(nodes) %>%
    dplyr::right_join(data)

  p <- ggraph::ggraph(
    arch, layout = "graphopt", charge = 0.1, mass = 30, niter = 100000
  ) +
    ggraph::geom_edge_link(
      ggplot2::aes_string(color = "trait"), width = 2, alpha = 0.6
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
