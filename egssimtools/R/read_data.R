#' Read data
#'
#' @param filename Path to the folder
#' @param variable What to read
#'
#' @export

read_data <- function(folder, variable) {

  read_binary(paste0(folder, "/", variable, ".dat"))

}
