#' Get locus chromosomal locations
#'
#' @param arch A genetic architecture
#'
#' @return A vector of chromosome indices for each locus
#'
#' @export

get_chromosomes <- function(arch) {

  library(tidyverse)
  arch$locations %>% map_int(~ min(which(.x <= arch$chromosomes)))

}
