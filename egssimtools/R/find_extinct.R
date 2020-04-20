#' Find extinct simulations
#'
#' @param root Path to the root directory of the simulations to check
#' @param pattern Pattern defining the simulation folders to look into
#' @param pb Whether to display a progress bar
#'
#' @export

# How many simulations did go extinct?
find_extinct <- function(root, pattern = "^sim_", pb = TRUE) {

  # Use progress bar?
  if (pb) {
    library(pbapply)
    myapply <- pbsapply
  } else myapply <- sapply

  # Get a list of folders in the directory
  folders <- list.files(root, pattern = pattern, full.names = TRUE)

  # Check each one for extinction
  extincts <- myapply(folders, is_extinct)

  # Return the names of the extinct folders
  if (any(extincts)) return (folders[extincts])
  return (NULL)

}
