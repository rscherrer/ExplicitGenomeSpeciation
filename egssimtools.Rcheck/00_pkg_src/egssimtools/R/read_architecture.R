#' Read architecture
#'
#' @param folder Path to the folder
#' @param filename Optional file name
#'
#' @export

read_architecture <- function(folder, filename = "architecture.txt") {

  read_archfile(paste0(folder, '/', filename))

}
