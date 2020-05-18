#' Find missing simulations
#'
#' @param simulations Either path to a root directory containing simulation folders, or a vector of simulation folders
#' @param pattern Pattern defining the simulation folders to look into
#' @param verbose Whether to display messages
#' @param pb Whether to display a progress bar
#' @param level Recursion level
#'
#' @export

# How many simulations are missing?
find_missing <- function(simulations, pattern = "^sim_", verbose = TRUE, pb = TRUE, level = 0) {

  library(pbapply)

  if (!verbose) pb <- FALSE else message("Looking for missing simulations...")
  if (pb) thissapply <- pbsapply else thissapply <- sapply
  if (level > 0) simulations <- fetch_dirs(simulations, pattern = pattern, level = level)

  # Check each one for missing data
  missings <- thissapply(simulations, is_missing)

  # Return the names of the extinct folders
  if (any(missings)) return (simulations[missings])
  return (NULL)

}
