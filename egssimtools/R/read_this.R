#' Wrapper around data readers
#'
#' This function works as an alias for the `read_*` family of functions
#'
#' @param root Path to a simulation folder
#' @param type Which `read_*` function to call, either of `data`, `sim`, `pop`,
#' `genome` or `network`
#' @param ... Parameters to be passed to the reading function
#'
#' @return A tibble
#'
#' @examples
#'
#' \dontrun{
#'
#' root <- "data/example_1"
#'
#' read_this(root, "sim", variables = "EI")
#' read_this(root, "data", variables = c("time", "EI"))
#'
#' }
#'
#' @seealso `read_data`, `read_sim`, `read_pop`, `read_genome`, `read_network`
#'
#' @export

read_this <- function(root, type = "data", ...) {
  if (type == "data") {
    read_data(root, ...)
  } else if (type == "sim") {
    read_sim(root, ...)
  } else if (type == "pop") {
    read_pop(root, ...)
  } else if (type == "genome") {
    read_genome(root, ...)
  } else if (type == "network") {
    read_network(root, ...)
  } else {
    stop("invalid value for type argument")
  }
}
