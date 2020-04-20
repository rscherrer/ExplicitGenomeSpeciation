#' Read ecotype trait means
#'
#' Means of the three traits within each ecotype for each time point
#'
#' @param filename Path to the folder
#'
#' @export

read_ecotype_means <- function(folder) {

  data <- read_data(folder, "ecotype_means")
  t <- read_time(folder)
  t <- rep(t, each = 6) # per time point
  lapply(split(data, t), split, c(1, 2, 3)) # per trait per ecotype

}
