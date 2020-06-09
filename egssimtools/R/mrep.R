#' Repeat many elements many times
#'
#' @param x A vector of things to repeat
#' @param n A vector of how many times to repeat each thing
#'
#' @return A vector of the same type as `x`

mrep <- function(x, n) {

  do.call("c", purrr::map2(x, n, ~ rep(.x, .y)))

}
