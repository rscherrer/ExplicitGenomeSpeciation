#' Read resource concentrations
#'
#' Equilibrium resource concentrations by resource by habitat by time point.
#'
#' @param filename Path to the folder
#'
#' @export

read_resources <- function(folder) {

  data <- read_data(folder, "resources")
  t <- read_time(folder)
  t <- rep(t, each = 4) # per time point
  lapply(split(data, t), split, c(1, 1, 2, 2)) # per habitat

}
