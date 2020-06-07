#' Get locus chromosomal locations
#'
#' @param arch A genetic architecture
#'
#' @return A vector of chromosome indices for each locus
#'
#' @export

get_chromosomes <- function(arch) {

   purrr::map_int(arch$locations, ~ min(which(.x <= arch$chromosomes)))

}
