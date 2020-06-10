#' Guess the number of rows per time point
#'
#' Function to find out the number of rows per time point in a given simulation
#'
#' @param folder Path to the folder
#' @param variable The variable to use to guess the number of rows
#'
#' @details The file `time.dat` must be present
#'
#' @return An integer
#'
#' @examples
#'
#' \dontrun{
#'
#' # Location of the simulation folder
#' root <- "data/example_1"
#'
#' guess_nrows(root, "EI")
#'
#' }
#'
#' @export

guess_nrows <- function(folder, variable) {

  time <- read_binary(paste0(folder, "/time.dat"))
  ntimes <- length(time)
  x <- read_binary(paste0(folder, "/", variable, ".dat"))
  length(x) / ntimes

}
