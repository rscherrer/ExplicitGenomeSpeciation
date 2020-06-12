#' Find extinct simulations
#'
#' @param sims A vector of simulation folders
#' @param ... Parameters for `is_extinct`
#'
#' @return A character vector of extinct simulations, or NULL if no simulation
#'   is extinct
#'
#' @examples
#'
#' \dontrun{
#'
#' # Location of the simulation folder
#' root <- "data"
#'
#' roots <- fetch_dirs(root, "example", level = 1)
#' find_extinct(roots)
#'
#' }
#'
#' @export

# How many simulations did go extinct?
find_extinct <- function(sims, ...) {

  # Check each one for extinction
  extincts <- purrr::map_lgl(sims, is_extinct, ...)

  # Return the names of the extinct folders
  if (any(extincts)) return (sims[extincts])
  return (NULL)

}

