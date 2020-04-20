#' Read individual data
#'
#' @param filename Path to the folder
#' @param variable What to read
#'
#' @export

read_individuals <- function(folder, variable) {

  t <- read_time(folder)
  n <- read_total_counts(folder)
  t_indiv <- mrep(t, n) # time at which each individual lived

  data <- read_data(folder, variable)
  split(data, t_indiv)

}
