#' Find missing simulations
#'
#' @param root Path to the root directory of the simulations to check
#' @param pb Whether to display a progress bar
#'
#' @export

# How many simulations are missing?
find_missing <- function(root, pb = TRUE) {

  # Use progress bar?
  if (pb) {
    library(pbapply)
    myapply <- pbsapply
  } else myapply <- sapply

  # Get a list of folders in the directory
  folders <- list.files(root, pattern = "^sim_", full.names = TRUE)

  # Check each one for missing data
  missings <- myapply(folders, is_missing)

  # Return the names of the extinct folders
  if (any(missings)) return (folders[missings])
  return (NULL)

}
