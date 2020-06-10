#' Get locus chromosomal locations
#'
#' @param arch A genetic architecture
#'
#' @return A vector of chromosome indices for each locus
#'
#' @examples
#'
#' \dontrun{
#'
#' # Location of the simulation folder
#' root <- "data/example_1"
#'
#' arch <- read_arch(root)
#' get_chromosomes(arch)
#'
#' }
#'
#' @export

get_chromosomes <- function(arch) {

   purrr::map_int(arch$locations, ~ min(which(.x <= arch$chromosomes)))

}
