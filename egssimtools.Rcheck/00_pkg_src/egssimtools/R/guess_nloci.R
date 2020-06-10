#' Guess the number of loci
#'
#' Function to find out the number of loci in a given simulation
#'
#' @param folder Path to the folder
#' @param variable The locus-specific variable to use to guess the number of loci
#'
#' @details The file `time.dat` must be present
#'
#' @return An integer
#'
#' @examples
#'
#' root <- system.file("extdata", "example_1", package = "egssimtools")
#' guess_nloci(root)
#'
#' @export

guess_nloci <- function(folder, variable = "genome_Fst") {

  time <- read_binary(paste0(folder, "/time.dat"))
  ntimes <- length(time)
  x <- read_binary(paste0(folder, "/", variable, ".dat"))
  length(x) / ntimes

}
