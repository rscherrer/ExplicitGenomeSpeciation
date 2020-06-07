#' Read architecture
#'
#' @param folder Path to the folder
#' @param filename Optional file name
#'
#' @return A list with all the fields of a genetic architecture
#'
#' @export

read_architecture <- function(folder, filename = "architecture.txt") {

  read_archfile(paste0(folder, '/', filename))

}
