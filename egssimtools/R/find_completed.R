#' Find completed simulations
#'
#' @param simulations Either path to a root directory containing simulation folders, or a vector of simulation folders
#' @param pattern Pattern defining the simulation folders to look into
#' @param verbose Whether to display messages
#' @param pb Whether to display a progress bar
#'
#' @export

find_completed <- function(simulations, pattern = "^sim_", verbose = TRUE, pb = TRUE) {

  library(pbapply)

  if (!verbose) pb <- FALSE else message("Looking for completed simulations...")

  # Get a list of folders in the directory
  if (length(simulations) == 1) simulations <- list.files(simulations, pattern = pattern, full.names = TRUE)

  status <- collect_status(simulations, verbose = verbose, pb = pb)
  if (length(which(status == "COMPLETED")) == 0) return (NULL)
  return (simulations[status == "COMPLETED"])

}

