#' Plot genomic trajectories
#'
#' Trajectories of loci through time along a certain variable.
#'
#' @param data A locus-wise data frame with columns "time", "locus" and `y` at least
#' @param y Name of the variable to plot
#' @param color Name of the variable to map to that aesthetic.
#' @param alpha Transparency value
#' @param by_trait Whether to facet by trait
#'
#' @return A ggplot
#'
#' @export

plot_genome_trajectories <- function(

  data,
  y,
  color = NULL,
  alpha = 0.5,
  by_trait = FALSE

) {

  # Set up
  p <- ggplot2::ggplot(data)

  # Prepare aesthetic mapping
  aes_string <- "x = time, y = get(y), group = locus"
  if (!is.null(color)) aes_string <- paste0(aes_string, ", color = get(color)")
  aes_string <- paste0("ggplot2::aes(", aes_string, ")")

  # Add lines
  p <- p + geom_line(mapping = eval(rlang::parse_expr(aes_string)), alpha = alpha)

  # Relabel axis
  p <- p + ylab(y)

  # Facettize if needed
  if (by_trait) p <- ggsim::facettize(p, cols = "trait", prepend = "Trait ")

  return(p)

}
