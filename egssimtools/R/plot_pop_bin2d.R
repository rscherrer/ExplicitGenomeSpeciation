#' Plot population distribution through time
#'
#' @param root Path to the simulation folder
#' @param y Variable to plot the density for
#' @param x Optional variable for the horizontal axis. Defaults to "time".
#' @param by Numbers of columns in which to split variable `x` (only if `x` is
#' not "time" and not `y`) and `y`. See `?read_data`.
#' @param j If variable `y` is splitted, which column to show?
#' @param k If variable `x` is splitted, which column to show?
#' @param t Optional time point to subset (makes sense if the heatmap is not
#' through time). Provide a negative number to automatically filter only the
#' last time point.
#' @param ... Parameters to be passed to `geom_path`
#'
#' @return A ggplot
#'
#' @examples
#'
#' root <- system.file("extdata", "example_1", package = "egssimtools")
#' plot_pop_bin2d(root, "individual_trait", by = 3, j = 1)
#'
#' # Plot the mating versus ecological trait
#' plot_pop_bin2d(
#'   root, "individual_trait", x = "individual_trait", by = 3,
#'   j = 2, k = 1
#' )
#'
#' @export

plot_pop_bin2d <- function(
  root,
  y,
  x = "time",
  by = 1,
  j = 1,
  k = 1,
  t = NULL,
  ...
) {

  to_read <- y
  if (x != "time" & x != y) to_read <- c(x, y)
  data <- read_indiv(root, to_read, by)

  xcol <- colnames(data)[grep(x, colnames(data))][k]
  ycol <- colnames(data)[grep(y, colnames(data))][j]

  if (!is.null(t)) {
    if (t < 0) t <- dplyr::last(data$time)
    data <- data %>% dplyr::filter(time == t)
  }

  ggplot2::ggplot(data, ggplot2::aes_string(x = xcol, y = ycol)) +
    ggplot2::geom_bin2d(...) +
    ggplot2::theme_bw()

}
