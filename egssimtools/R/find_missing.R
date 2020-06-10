#' Find missing simulations
#'
#' @param sims Either path to a root directory containing simulation folders, or
#'   a vector of simulation folders
#' @param pattern Pattern defining the simulation folders to look into
#' @param verbose Whether to display messages
#' @param level Recursion level
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
#' find_missing(root, pattern = "example_", level = 1)
#'
#' }
#'
#' @export

# How many simulations are missing?
find_missing <- function(
  sims,
  pattern = "^sim_",
  verbose = TRUE,
  level = 0
) {

  if (verbose) message("Looking for missing simulations...")
  if (level > 0) sims <- fetch_dirs(sims, pattern = pattern, level = level)

  # Check each one for missing data
  missings <- sapply(sims, is_missing)

  # Return the names of the extinct folders
  if (any(missings)) return (sims[missings])
  return (NULL)

}
