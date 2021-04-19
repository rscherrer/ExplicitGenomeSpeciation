#' Collect SLURM exit status
#'
#' @param simulations Either path to a root directory containing simulation folders, or a vector of simulation folders
#' @param pattern Pattern defining the simulation folders to look into
#' @param verbose Whether to display messages
#' @param level Recursion level
#'
#' @return A character vector
#'
#' @export

collect_status <- function(
  simulations,
  pattern = "^sim_",
  verbose = TRUE,
  level = 0
) {

  if (verbose) message("Collecting SLURM exit status...")
  if (level > 0) simulations <- fetch_dirs(
    simulations, pattern = pattern, level = level
  )
  return (sapply(simulations, get_status))

}
