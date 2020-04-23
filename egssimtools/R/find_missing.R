#' Find missing simulations
#'
#' @param simulations Either path to a root directory containing simulation folders, or a vector of simulation folders
#' @param pattern Pattern defining the simulation folders to look into
#' @param pb Whether to display a progress bar
#'
#' @export

# How many simulations are missing?
find_missing <- function(simulations, pattern = "^sim_", pb = TRUE) {

  library(pbapply)

  # Use progress bar?
  if (pb) thissapply <- pbsapply else thissapply <- sapply

  # Get a list of folders in the directory
  if (length(simulations) == 1) simulations <- list.files(simulations, pattern = pattern, full.names = TRUE)

  # Check each one for missing data
  missings <- thissapply(simulations, is_missing)

  # Return the names of the extinct folders
  if (any(missings)) return (simulations[missings])
  return (NULL)

}
