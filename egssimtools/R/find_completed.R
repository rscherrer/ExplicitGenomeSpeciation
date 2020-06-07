#' Find completed simulations
#'
#' @param sims Either path to a root directory containing simulation folders, or a vector of simulation folders
#' @param pattern Pattern defining the simulation folders to look into
#' @param verbose Whether to display messages
#' @param level Recursion level
#'
#' @return A character vector with the paths to completed simulations
#'
#' @export

find_completed <- function(
  sims,
  pattern = "^sim_",
  verbose = TRUE,
  level = 0
) {

  if (verbose) message("Looking for completed simulations...")
  if (level > 0) sims <- fetch_dirs(
    sims, pattern = pattern, level = level
  )

  status <- collect_status(sims, verbose = verbose)
  if (length(which(status == "COMPLETED")) == 0) return (NULL)
  return (sims[status == "COMPLETED"])

}

