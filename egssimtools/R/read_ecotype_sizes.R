#' Read ecotype population size
#'
#' Population size through time within each ecotype
#'
#' @param filename Path to the folder
#'
#' @export

read_ecotype_sizes <- function(folder) {

  data <- read_data(folder, variable = "ecotype_size")
  split(data, rep(seq_len(length(data) / 2), each = 2))

}
