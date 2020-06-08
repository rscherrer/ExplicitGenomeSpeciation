#' Diagnostic genome scan
#'
#' Quick-plots a genome scan of a given variable at a given time point
#'
#' @param root Path to the simulation folder
#' @param y Name of the variable to plot
#' @param x Optional x-axis, defaults to locus identity but could also be e.g. "location" for locus position along the genome
#' @param t Time point at which to show the scan. If unspecified, last time point.
#'
#' @return A ggplot
#'
#' @export

plot_genome_scan <- function(root, y, x = "locus", t = NULL) {

  data <- read_loci(root, y, architecture = TRUE, archfile = archfile)

  if (is.null(t)) t <- dplyr::last(data$time)

  colors <- c("forestgreen", "goldenrod", "lightgrey")

  data <- data %>% dplyr::filter(time == t)
  ggplot2::ggplot(
    data,
    ggplot2::aes(x = get(x), y = get(y), color = factor(trait))
  ) +
    ggplot2::geom_point() +
    ggplot2::theme_bw() +
    ggplot2::labs(x = x, y = y, color = "Trait") +
    ggplot2::scale_color_manual(values = colors)

}
