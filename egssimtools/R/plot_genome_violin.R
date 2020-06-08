#' Diagnostic violin plot
#'
#' Quick-plots a variable across the genome, possibly split between categories
#'
#' @param root Path to the simulation folder
#' @param y Name of the variable to plot
#' @param x Optional grouping variable
#' @param t Time point at which to show the scan. If unspecified, last time point.
#' @param archfile Optional architecture file name
#'
#' @return A ggplot
#'
#' @export

plot_genome_violin <- function(
  root,
  y,
  x = NULL,
  t = NULL,
  archfile = "architecture.txt"
) {

  data <- read_loci(root, y, architecture = TRUE, archfile = archfile)

  if (is.null(t)) t <- dplyr::last(data$time)

  data <- data %>% dplyr::filter(time == t)

  colors <- c("forestgreen", "goldenrod", "lightgrey")

  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = factor(get(x)), y = get(y), color = factor(trait))
  ) +
    ggplot2::geom_violin() +
    ggplot2::theme_bw() +
    ggplot2::labs(x = x, y = y) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(x = "Trait", y = y, color = "Trait")

  return(p)

}
