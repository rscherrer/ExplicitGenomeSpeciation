#' Plot a simulation-wise variable
#'
#' @param root Path to the simulation folder
#' @param y Variable to plot
#' @param x Optional variable for the horizontal axis. Defaults to "time".
#' @param by Numbers of columns in which to split variable `x` (only if `x` is not "time") and `y`. See `?read_data`.
#' @param j If variable `y` is splitted, which column to show?
#' @param k If variable `x` is splitted, which column to show?
#' @param ... Parameters to be passed to `geom_path`
#'
#' @return A ggplot
#'
#' @examples
#'
#' root <- system.file("extdata", "example_1", package = "egssimtools")
#' plot_sim_line(root, "EI")
#' plot_sim_line(root, "RI", x = "EI")
#'
#' @export

plot_sim_line <- function(root, y, x = "time", by = NULL, j = 1, k = 1, ...) {

  to_read <- y
  if (x != "time" & x != y) to_read <- c(x, y)
  if (is.null(by)) by <- rep(1, length(to_read))
  data <- read_sim(root, to_read, by)

  xcol <- colnames(data)[grep(x, colnames(data))][k]
  ycol <- colnames(data)[grep(y, colnames(data))][j]

  ggplot2::ggplot(data, ggplot2::aes(x = get(xcol), y = get(ycol))) +
    ggplot2::geom_path(...) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = xcol, y = ycol)

}
