#' Find extant simulations
#'
#' Find simulations that did not crash, that successfully completed and that did
#' not go extinct
#'
#' @param roots Vector of simulation folders
#' @param ... Parameters for `is_extant`
#'
#' @return A character vector
#'
#' @examples
#'
#' \dontrun{
#'
#' roots <- fetch_dirs("data", "example", 1)
#' find_extant(roots)
#'
#' }
#'
#' @export

find_extant <- function(roots, ...) {

  extant <- purrr::map_lgl(roots, is_extant, ...)
  return(roots[extant])

}
