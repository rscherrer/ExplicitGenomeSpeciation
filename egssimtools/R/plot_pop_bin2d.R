#' Plot population distribution through time
#'
#' @param root Path to the simulation folder
#' @param y Variable to plot the density for
#' @param by How many columns to split the variable into? See `?read_data`.
#' @param j If the variable is splitted, which column to show?
#' @param bins Number of bins for `geom_bin2d`
#'
#' @return A ggplot
#'
#' @examples
#'
#' root <- system.file("extdata", "example_1", package = "egssimtools")
#' plot_pop_bin2d(root, "individual_trait", by = 3, j = 1)
#'
#' @export

plot_pop_bin2d <- function(root, y, by = 1, j = 1, bins = 50) {

  data <- read_indiv(root, y, by)
  data <- data[, c(1, j + 1)]

  if (is.null(t)) t <- dplyr::last(data$time)

  ggplot2::ggplot(data, ggplot2::aes(x = time, y = get(colnames(data)[2]))) +
    ggplot2::geom_bin2d(bins = bins) +
    ggplot2::theme_bw() +
    ggplot2::ylab(colnames(data)[2])

}
