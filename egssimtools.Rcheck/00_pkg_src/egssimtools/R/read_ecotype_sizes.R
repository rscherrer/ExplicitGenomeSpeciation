#' Read ecotype population size
#'
#' Population size through time within each ecotype
#'
#' @param filename Path to the folder
#'
#' @export

read_ecotype_sizes <- function(folder) {

  data <- read_data(folder, variable = "ecotype_size")
  t <- read_time_if(folder, ntimes = 2)
  split(data, rep(t, each = 2))

}
