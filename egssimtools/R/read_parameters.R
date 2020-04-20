#' Read parameters
#'
#' @param folder Path to the folder
#' @param parnames Optional vector of parameter names
#' @param filename Optional file name
#'
#' @export

read_parameters <- function(folder, parnames = NULL, filename = "paramlog.txt") {

  read_paramfile(paste0(folder, '/', filename), parnames)

}
