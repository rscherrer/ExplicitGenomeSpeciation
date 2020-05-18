#' Find extinct simulations
#'
#' @param simulations Either path to a root directory containing simulation folders, or a vector of simulation folders
#' @param pattern Pattern defining the simulation folders to look into
#' @param verbose Whether to display messages
#' @param pb Whether to display a progress bar
#' @param level Recursion level
#'
#' @export

# How many simulations did go extinct?
find_extinct <- function(simulations, pattern = "^sim_", verbose = TRUE, pb = TRUE, level = 0) {

  library(pbapply)

  if (!verbose) pb <- FALSE else message("Looking for extinct simulations...")
  if (pb) thissapply <- pbsapply else thissapply <- sapply
  if (level > 0) simulations <- fetch_dirs(simulations, pattern = pattern, level = level)

  # Check each one for extinction
  extincts <- thissapply(simulations, is_extinct)

  # Return the names of the extinct folders
  if (any(extincts)) return (simulations[extincts])
  return (NULL)

}

