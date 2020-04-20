#' Read trait means
#'
#' Means of the three traits for each time point
#'
#' @param filename Path to the folder
#'
#' @export

read_means <- function(folder) {

  data <- read_data(folder, "means")
  t <- read_time(folder)
  t <- rep(t, each = 3) # per time point
  split(data, t)

}
