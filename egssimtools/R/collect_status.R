#' Collect SLURM exit status
#'
#' @param simulations Either path to a root directory containing simulation folders, or a vector of simulation folders
#' @param pattern Pattern defining the simulation folders to look into
#' @param verbose Whether to display messages
#' @param pb Whether to display a progress bar
#'
#' @export

collect_status <- function(simulations, pattern = "^sim_", verbose = TRUE, pb = TRUE) {

  library(pbapply)

  if (!verbose) pb <- FALSE else message("Collecting SLURM exit status...")

  # Use progress bar?
  if (pb) thissapply <- pbsapply else thissapply <- sapply

  # Get a list of folders in the directory
  if (length(simulations) == 1) simulations <- list.files(simulations, pattern = pattern, full.names = TRUE)

  # Check each one for extinction
  return (thissapply(simulations, get_status))

}
