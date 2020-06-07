#' Diagnose density of a given variable across the population
#'
#' @param root Path to the simulation folder
#' @param y Variable to plot the density for
#' @param by How many columns to split the variable into? See `?read_data`.
#' @param j If the variable is splitted, which column to show?
#' @param t Optional time point to filter. Defaults to the last time point
#' @param fill Optional fill color
#'
#' @return A ggplot
#'
#' @export

plot_population_density <- function(root, y, by = 1, j = 1, t = NULL, fill = "lightgrey") {

  library(ggplot2)

  data <- read_data(root, c("time", y), by = c(1, by), dupl = list("population_size", 1))
  data <- data[, c(1, j + 1)]

  if (is.null(t)) t <- last(data$time)

  data <- data %>% filter(time == t)

  ggplot(data, aes(x = get(colnames(data)[2]))) +
    geom_density(fill = fill) +
    theme_bw() +
    xlab(colnames(data)[2])

}
