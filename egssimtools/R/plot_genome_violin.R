#' Diagnostic violin plot
#'
#' Quick-plots a variable across the genome, possibly split between categories
#'
#' @param root Path to the simulation folder
#' @param y Name of the variable to plot
#' @param x Optional grouping variable
#' @param t Time point at which to show the scan. If unspecified, last time point.
#' @param is_color Whether to color by `x`
#' @param archfile Optional architecture file name
#'
#' @return A ggplot
#'
#' @examples
#'
#' root <- system.file("extdata", "example_1", package = "egssimtools")
#' plot_genome_violin(root, "genome_Fst")
#' plot_genome_violin(root, "genome_Fst", x = "trait", is_color = TRUE)
#'
#' @export

plot_genome_violin <- function(
  root,
  y,
  x = NULL,
  t = NULL,
  is_color = FALSE,
  archfile = "architecture.txt"
) {

  data <- read_loci(root, y, architecture = TRUE, archfile = archfile)

  if (is.null(t)) t <- dplyr::last(data$time)

  data <- data %>% dplyr::filter(time == t)

  colors <- c("forestgreen", "goldenrod", "lightgrey")

  if (is.null(x)) x <- "x"

  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = factor(get(x)), y = get(y))
  )
  if (is_color) p <- p + ggplot2::aes(color = factor(get(x)))
  p <- p +
    ggplot2::geom_violin() +
    ggplot2::theme_bw() +
    ggplot2::labs(x = x, y = y) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(x = x, y = y, color = x)
  if (x == "x") p <- p +
    ggplot2::xlab(NULL) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank())

  return(p)

}
