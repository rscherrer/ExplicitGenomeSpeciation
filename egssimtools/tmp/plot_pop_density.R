#' Plot population distribution
#'
#' @param root Path to the simulation folder
#' @param y Variable to plot the density for
#' @param by How many columns to split the variable into? See `?read_data`.
#' @param j If the variable is splitted, which column to show?
#' @param plot_type See `?ggsim::ggdensityplot`
#' @param t Optional time point to filter. Defaults to the last time point
#' @param fill Optional fill color
#' @param ... Parameters to be passed to `ggsim::ggdensityplot`
#'
#' @return A ggplot
#'
#' @examples
#'
#' # Plot the distribution of the ecological trait
#' root <- system.file("extdata", "example_1", package = "egssimtools")
#' plot_pop_density(root, "individual_trait", by = 3, j = 1)
#'
#' @export

plot_pop_density <- function(
  root,
  y,
  by = 1,
  j = 1,
  plot_type = "density",
  t = NULL,
  ...
) {

  data <- read_indiv(root, y, by)
  data <- data[, c(1, j + 1)]
  if (is.null(t)) t <- dplyr::last(data$time)
  data <- data %>% dplyr::filter(time == t)
  ggsim::ggdensityplot(data, colnames(data)[2], plot_type, ...)

}
