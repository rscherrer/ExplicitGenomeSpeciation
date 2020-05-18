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

dplot_genome_scan <- function(root, y, x = "locus", t = NULL) {

  library(ggplot2)

  time <- read_data(root, "time")
  X <- read_data(root, y)
  nloci <- nrow(X) / nrow(time)

  data <- read_data(root, c("time", y), dupl = c(nloci, 1), architecture = TRUE)

  if (is.null(t)) t <- last(data$time)

  data <- data %>% filter(time == t)
  ggplot(data, aes(x = get(x), y = get(y))) +
    geom_point() +
    theme_bw() +
    labs(x = x, y = y)

}
