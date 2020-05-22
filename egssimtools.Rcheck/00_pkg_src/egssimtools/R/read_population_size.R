#' Read total population size
#'
#' @param filename Path to the folder
#'
#' @export

read_population_size <- function(folder) {

  read_data(folder, variable = "population_size")

}
