#' Guess the number of loci
#'
#' Function to find out the number of loci in a given simulation
#'
#' @param folder Path to the folder
#' @param variable The locus-specific variable to use to guess the number of
#' loci
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
#' guess_nloci(root)
#'
#' }
#'
#' @export

guess_nloci <- function(folder, variable = "genome_Fst") {

  guess_nrows(folder, variable)

}
