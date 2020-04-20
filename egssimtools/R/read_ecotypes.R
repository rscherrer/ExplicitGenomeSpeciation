#' Read individual ecotypes
#'
#' @param filename Path to the folder
#'
#' @export

read_ecotypes <- function(folder) {

  read_individuals(folder, "individual_ecotype", nvalues = 1)

}
