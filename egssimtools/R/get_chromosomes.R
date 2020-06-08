#' Get locus chromosomal locations
#'
#' @param arch A genetic architecture
#'
#' @return A vector of chromosome indices for each locus
#'
#' @examples
#'
#' root <- system.file("extdata", "example_1", package = "egssimtools")
#' arch <- read_architecture(root)
#' get_chromosomes(arch)
#'
#' @export

get_chromosomes <- function(arch) {

   purrr::map_int(arch$locations, ~ min(which(.x <= arch$chromosomes)))

}
