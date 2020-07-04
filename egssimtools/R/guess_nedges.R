#' Guess the number of edges
#'
#' Function to find out the number of edges in a given simulation
#'
#' @param folder Path to the folder
#' @param variable The edge-specific variable to use to guess the number of
#' edges
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
#' guess_nedges(root)
#'
#' }
#'
#' @export

guess_nedges <- function(folder, variable = "network_corfreq") {

  guess_nrows(folder, variable)

}
