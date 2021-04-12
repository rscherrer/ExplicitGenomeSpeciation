#' Find the partners of a locus in a network
#'
#' @param locus The focal locus
#' @param arch The network-side of the architecture (second list element
#' returned by `read_arch_network`, with `as_list` set to TRUE)
#'
#' @return A vector of integers
#'
#' @export

find_partners <- function(locus, arch) {

  partners <- arch %>% dplyr::filter(from == locus | to == locus)
  partners <- c(partners$from, partners$to)
  return(partners[partners != locus])

}
