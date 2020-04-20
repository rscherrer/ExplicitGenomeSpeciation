#' Read individual habitats
#'
#' @param filename Path to the folder
#'
#' @export

read_habitats <- function(folder) {

  read_individuals(folder, "individual_habitat", nvalues = 1)

}
