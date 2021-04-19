#' Plot genome scan
#'
#' Plots a genome scan of a given variable at a given time point
#'
#' @param root Path to the simulation folder
#' @param y Name of the variable to plot
#' @param x Optional x-axis, defaults to locus identity but could also be e.g.
#' "location" for locus position along the genome
#' @param t Time point at which to show the scan. If unspecified, last time
#' point.
#' @param mapping Named list of extra mappings to be passed to
#' `ggplot2::aes_string`
#' @param archfile Optional architecture file name
#' @param ... Extra parameters for `ggplot2::geom_point`
#'
#' @return A ggplot
#'
#' @examples
#'
#' root <- system.file("extdata", "example_1", package = "egssimtools")
#' plot_genome_scan(root, "genome_Fst", mapping = list(color = "trait"))
#'
#' @export

plot_genome_scan <- function(
  root,
  y,
  x = "locus",
  t = NULL,
  mapping = NULL,
  archfile = "architecture.txt",
  ...
) {

  data <- read_loci(root, y, architecture = TRUE, archfile = archfile)

  if (is.null(t)) t <- dplyr::last(data$time)

  data <- data %>% dplyr::filter(time == t)
  p <- ggplot2::ggplot(
    data,
    ggplot2::aes_string(x = x, y = y)
  )
  if (!is.null(mapping)) {
    this_aes <- ggplot2::aes_string
    this_mapping <- do.call("this_aes", mapping)
    p <- p + this_mapping
  }
  p <- p +
    ggplot2::geom_point(...) +
    ggplot2::theme_bw()
  p

}
