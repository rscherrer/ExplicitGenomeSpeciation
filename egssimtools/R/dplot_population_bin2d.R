#' Diagnose density of a given variable across the population through time
#'
#' @param root Path to the simulation folder
#' @param y Variable to plot the density for
#' @param by How many columns to split the variable into? See `?read_data`.
#' @param j If the variable is splitted, which column to show?
#' @param bins Number of bins for `geom_bin2d`
#'
#' @return A ggplot
#'
#' @export

dplot_population_bin2d <- function(root, y, by = 1, j = 1, bins = 50) {

  library(ggplot2)

  data <- read_data(root, c("time", y), by = c(1, by), dupl = list("population_size", 1))
  data <- data[, c(1, j + 1)]

  if (is.null(t)) t <- last(data$time)

  ggplot(data, aes(x = time, y = get(colnames(data)[2]))) +
    geom_bin2d(bins = bins) +
    theme_bw() +
    ylab(colnames(data)[2])

}
