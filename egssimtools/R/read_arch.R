#' Read architecture
#'
#' @param folder Path to the folder
#' @param filename Optional file name
#'
#' @return A list with all the fields of a genetic architecture
#'
#' @examples
#'
#' \dontrun{
#'
#' # Location of the simulation folder
#' root <- "data/example_1"
#'
#' read_arch(root)
#'
#' }
#'
#' @export

read_arch <- function(folder, filename = "architecture.txt") {

  read_arch_file(paste0(folder, '/', filename))

}
