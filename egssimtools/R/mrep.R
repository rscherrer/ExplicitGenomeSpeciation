#' Repeat multiple elements multiple times
#'
#' @param x Vector of what to repeat
#' @param n Vector of how many times
#'
#' @export

mrep <- function(x, n) {

  if (all(n == n[1])) return (rep(x, each = n[1]))
  do.call("c", mapply(rep, x, n))

}
