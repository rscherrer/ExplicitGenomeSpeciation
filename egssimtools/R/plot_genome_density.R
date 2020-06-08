#' Plot genomic distribution
#'
#' @param root Path to the simulation folder
#' @param y Name of the variable to plot
#' @param plot_type See `?ggsim::ggdensityplot`
#' @param t Optional time point. Last time point if unspecified.
#' @param archfile Optional architecture file name
#' @param ... Parameters to be passed to `ggsim::ggdensityplot`
#'
#' @return A ggplot
#'
#' @examples
#'
#' root <- system.file("extdata", "example_1", package = "egssimtools")
#' plot_genome_density(root, "genome_Fst", x = "trait")
#'
#' @export

plot_genome_density <- function(
  root,
  y,
  plot_type = "violin",
  t = NULL,
  archfile = "architecture.txt",
  ...
) {

  data <- read_loci(root, y, architecture = TRUE, archfile = archfile)
  if (is.null(t)) t <- dplyr::last(data$time)
  data <- data %>% dplyr::filter(time == t)
  ggsim::ggdensityplot(data, y, plot_type, ...)

}
