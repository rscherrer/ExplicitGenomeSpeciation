#' Collect SLURM exit status
#'
#' @param simulations Either path to a root directory containing simulation folders, or a vector of simulation folders
#' @param pattern Pattern defining the simulation folders to look into
#' @param verbose Whether to display messages
#' @param pb Whether to display a progress bar
#' @param level Recursion level
#'
#' @export

collect_status <- function(simulations, pattern = "^sim_", verbose = TRUE, pb = TRUE, level = 0) {

  library(pbapply)

  if (!verbose) pb <- FALSE else message("Collecting SLURM exit status...")
  if (pb) thissapply <- pbsapply else thissapply <- sapply
  if (level > 0) simulations <- fetch_dirs(simulations, pattern = pattern, level = level)

  # Check each one for extinction
  return (thissapply(simulations, get_status))

}
