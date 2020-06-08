#' Genome heatmap through time
#'
#' Quick-plots a given variable across the genome and through time
#'
#' @param root Path to the simulation folder
#' @param y Name of the variable to plot.
#' @param archfile The architecture file to look for
#'
#' @return A ggplot
#'
#' @examples
#'
#' root <- system.file("extdata", "example_1", package = "egssimtools")
#' plot_genome_heatmap(root, "genome_Fst")
#'
#' @export

plot_genome_heatmap <- function(root, y, archfile = "architecture.txt") {

  data <- read_loci(root, y, architecture = TRUE, archfile = archfile)

  ggplot2::ggplot(data, ggplot2::aes(x = time, y = locus, fill = get(y))) +
    ggplot2::geom_tile() +
    ggplot2::theme_bw() +
    ggplot2::labs(fill = y) +
    ggplot2::labs(x = "time", y = "locus", fill = y)

}
