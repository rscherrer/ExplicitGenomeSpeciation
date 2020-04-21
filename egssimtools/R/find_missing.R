#' Find missing simulations
#'
#' @param root Path to the root directory of the simulations to check
#' @param pattern Pattern defining the simulation folders to look into
#' @param pb Whether to display a progress bar
#'
#' @export

# How many simulations are missing?
find_missing <- function(root, pattern = "^sim_", pb = TRUE) {

  library(pbapply)

  # Use progress bar?
  if (pb) thissapply <- pbsapply else thissapply <- sapply

  # Get a list of folders in the directory
  folders <- list.files(root, pattern = pattern, full.names = TRUE)

  # Check each one for missing data
  missings <- thissapply(folders, is_missing)

  # Return the names of the extinct folders
  if (any(missings)) return (folders[missings])
  return (NULL)

}
