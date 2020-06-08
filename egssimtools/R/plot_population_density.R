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
#' @examples
#'
#' # Plot the distribution of the ecological trait
#' root <- system.file("extdata", "example_1", package = "egssimtools")
#' plot_population_density(root, "individual_trait", by = 3, j = 1)
#'
#' @export

plot_population_density <- function(
  root,
  y,
  by = 1,
  j = 1,
  t = NULL,
  fill = "lightgrey"
) {

  data <- read_indiv(root, y, by)
  data <- data[, c(1, j + 1)]

  if (is.null(t)) t <- dplyr::last(data$time)

  data <- data %>% dplyr::filter(time == t)

  ggplot2::ggplot(data, ggplot2::aes(x = get(colnames(data)[2]))) +
    ggplot2::geom_density(fill = fill) +
    ggplot2::theme_bw() +
    ggplot2::xlab(colnames(data)[2])

}
