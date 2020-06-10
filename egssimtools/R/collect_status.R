#' Collect SLURM exit status
#'
#' @param simulations Either path to a root directory containing simulation
#' folders, or a vector of simulation folders
#' @param pattern Pattern defining the simulation folders to look into
#' @param verbose Whether to display messages
#' @param level Recursion level
#'
#' @return A character vector
#'
#' @examples
#'
#' \dontrun{
#'
#' # Location of the simulation folder
#' root <- "data"
#'
#' # Should return not found with the example data because these simulations
#' # were run locally, not on SLURM
#' collect_status(root, pattern = "example", level = 1)
#'
#' }
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
