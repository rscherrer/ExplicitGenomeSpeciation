#' Find extant simulations
#'
#' @param simulations Either path to a root directory containing simulation folders, or a vector of simulation folders
#' @param pattern Pattern defining the simulation folders to look into
#' @param verbose Whether to display messages
#' @param pb Whether to display a progress bar
#'
#' @return A vector of names of the simulation folders that completed
#'
#' @export

find_extant <- function(simulations, pattern = "^sim_", verbose = TRUE, pb = TRUE) {

  library(pbapply)

  if (!verbose) pb <- FALSE else message("Looking for extant simulations...")
  if (pb) thislapply <- pblapply else thislapply <- lapply

  if (length(simulations) == 1) simulations <- list.files(simulations, pattern = pattern, full.names = TRUE)

  # Look for missing and extinct simulations
  missings <- find_missing(simulations, verbose = verbose, pb = pb)
  completed <- find_completed(simulations, verbose = verbose, pb = pb)
  extincts <- find_extinct(simulations, verbose = verbose, pb = pb)

  # Identify extant simulations
  if (!is.null(missings)) simulations <- simulations[!simulations %in% missings]
  if (!is.null(completed)) simulations <- simulations[simulations %in% completed]
  if (!is.null(extincts)) simulations <- simulations[!simulations %in% extincts]

  return (simulations)

}
