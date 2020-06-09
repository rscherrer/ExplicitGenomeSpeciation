#' Find completed simulations
#'
#' @param sims Either path to a root directory containing simulation folders, or a vector of simulation folders
#' @param pattern Pattern defining the simulation folders to look into
#' @param verbose Whether to display messages
#' @param level Recursion level
#'
#' @details Finds completed simulations by checking the SLURM exit status.
#' Applicable only to simulations run on a server with SLURM.
#'
#' @return A character vector with the paths to completed simulations
#'
#' @examples
#'
#' # Location of the simulation folder
#' root <- "egsimtools/data"
#'
#' find_completed(root, pattern = "example_", level = 1)
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

