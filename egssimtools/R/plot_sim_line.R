#' Plot a simulation-wise variable
#'
#' @param root Path to the simulation folder
#' @param y Variable to track
#' @param x Optional variable for the horizontal axis. Defaults to "time".
#' @param j If the variable is splitted, which column to show?
#' @param bins Number of bins for `geom_bin2d`
#' @param color Optional color of the line
#'
#' @return A ggplot
#'
#' @examples
#'
#' root <- system.file("extdata", "example_1", package = "egssimtools")
#' plot_sim_line(root, "EI")
#'
#' @export

plot_sim_line <- function(root, y, x = "time", by = 1, j = 1, color = "black") {

  data <- read_sim(root, y, by)
  data <- data[, c(1, j + 1)]

  if (is.null(t)) t <- dplyr::last(data$time)

  ggplot2::ggplot(data, ggplot2::aes(x = get(x), y = get(colnames(data)[2]))) +
    ggplot2::geom_path(color = color) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = x, y = colnames(data)[2])

}
