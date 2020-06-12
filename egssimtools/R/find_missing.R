#' Find missing simulations
#'
#' @param sims A vector of simulation folders
#'
#' @details Uses the absence of `.dat` output files as a clue for missing
#' simulations
#'
#' @return A character vector of missing simulations, or NULL if no simulation
#'   is missing
#'
#' @examples
#'
#' \dontrun{
#'
#' # Location of the simulation folder
#' root <- "data"
#'
#' roots <- fetch_dirs(root, "example_", level = 1)
#' find_missing(roots)
#'
#' }
#'
#' @export

# How many simulations are missing?
find_missing <- function(sims) {

  # Check each one for missing data
  missings <- purrr::map_lgl(sims, is_missing)

  # Return the names of the extinct folders
  if (any(missings)) return (sims[missings])
  return (NULL)

}
