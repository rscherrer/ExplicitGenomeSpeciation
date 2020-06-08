#' Plot gene network
#'
#' @param root Path to the simulation folder
#' @param y Name of the variable to map as node size aesthetics
#' @param t What time point? Defaults to last time point
#' @param archfile Optional architecture file name
#'
#' @return A ggplot
#'
#' @export

plot_network <- function(root, y, t = NULL, archfile = "architecture.txt") {

  data <- read_loci(root, y, architecture = TRUE, archfile = archfile)

  if (is.null(t)) t <- dplyr::last(data$time)
  data <- data %>% dplyr::filter(time == t)

  arch <- read_network_architecture(root)
  arch <- arch %>%
    tidygraph::activate(nodes) %>%
    dplyr::right_join(data) %>%
    tidygraph::activate(nodes) %>%
    dplyr::mutate(trait = factor(trait)) %>%
    tidygraph::activate(edges) %>%
    dplyr::mutate(trait = factor(trait))

  colors <- c("forestgreen", "goldenrod", "grey")
  p <- ggraph::ggraph(
    arch, layout = "graphopt", charge = 0.1, mass = 30, niter = 100000
  ) +
    ggraph::geom_edge_link(
      ggplot2::aes(color = factor(trait)), width = 2, alpha = 0.6
    ) +
    ggraph::geom_node_point(
      ggplot2::aes_string(fill = "trait", size = y), shape = 21
    ) +
    ggraph::scale_edge_color_manual(values = colors) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_size_continuous(range = c(2, 7)) +
    ggplot2::scale_alpha(range = c(0.6, 1)) +
    ggplot2::labs(size = y, fill = "Trait") +
    ggplot2::theme_void() +
    ggplot2::guides(edge_color = FALSE)

  p

}
