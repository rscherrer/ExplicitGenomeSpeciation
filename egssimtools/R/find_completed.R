#' Find completed simulations
#'
#' @param simulations Either path to a root directory containing simulation folders, or a vector of simulation folders
#' @param pattern Pattern defining the simulation folders to look into
#' @param verbose Whether to display messages
#' @param pb Whether to display a progress bar
#' @param level Recursion level
#'
#' @export

find_completed <- function(simulations, pattern = "^sim_", verbose = TRUE, pb = TRUE, level = 0) {

  library(pbapply)

  if (!verbose) pb <- FALSE else message("Looking for completed simulations...")
  if (level > 0) simulations <- fetch_dirs(simulations, pattern = pattern, level = level)

  status <- collect_status(simulations, verbose = verbose, pb = pb)
  if (length(which(status == "COMPLETED")) == 0) return (NULL)
  return (simulations[status == "COMPLETED"])

}

