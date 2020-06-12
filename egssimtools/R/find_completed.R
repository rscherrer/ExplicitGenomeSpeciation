#' Find completed simulations
#'
#' @param sims A vector of simulation folders
#'
#' @details Finds completed simulations by checking the SLURM exit status.
#' Applicable only to simulations run on a server with SLURM.
#'
#' @return A character vector with the paths to completed simulations
#'
#' @examples
#'
#' \dontrun{
#'
#' # Location of the simulation folder
#' root <- "data"
#'
#' roots <- fetch_dirs(roots, pattern = "example_", level = 1)
#' find_completed(roots)
#'
#' }
#'
#' @export

find_completed <- function(sims) {

  status <- collect_status(sims)
  if (length(which(status == "COMPLETED")) == 0) return (NULL)
  return (sims[status == "COMPLETED"])

}

