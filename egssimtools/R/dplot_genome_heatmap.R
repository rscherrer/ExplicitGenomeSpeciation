#' Diagnostic genome heatmap through time
#'
#' Quick-plots a given variable across the genome and through time
#'
#' @param root Path to the simulation folder
#' @param y Name of the variable to plot.
#' @param archfile The architecture file to look for. Defaults to "architecture.txt"
#'
#' @return A ggplot
#'
#' @export

dplot_genome_heatmap <- function(root, y, archfile = "architecture.txt") {

  library(ggplot2)

  time <- read_data(root, "time")
  X <- read_data(root, y)
  nloci <- nrow(X) / nrow(time)

  data <- read_data(
    root, c("time", y), dupl = c(nloci, 1),
    architecture = TRUE,
    archfile = archfile
  )

  ggplot(data, aes(x = time, y = locus, fill = get(y))) +
    geom_tile() +
    theme_bw() +
    labs(fill = y) +
    scale_fill_gradient(low = "black", high = "yellow") +
    labs(x = "Time", y = "Locus", fill = y)

}
