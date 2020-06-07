#' Diagnose a simulation variable
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
#' @export

plot_simulation_line <- function(root, y, x = "time", by = 1, j = 1, color = "black") {

  library(ggplot2)

  data <- read_data(root, c(x, y), by = c(1, by))
  data <- data[, c(1, j + 1)]

  if (is.null(t)) t <- last(data$time)

  ggplot(data, aes(x = get(x), y = get(colnames(data)[2]))) +
    geom_path(color = color) +
    theme_bw() +
    labs(x = x, y = colnames(data)[2])

}
