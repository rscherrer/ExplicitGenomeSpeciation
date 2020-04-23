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

  if (!verbose) pb <- FALSE
  if (pb) thislapply <- pblapply else thislapply <- lapply

  # Look for missing and extinct simulations
  if (verbose) message("Looking for missing simulations...")
  missings <- find_missing(simulations, pattern = pattern, pb = pb)
  if (verbose) message("Looking for extinct simulations...")
  extincts <- find_extinct(simulations, pattern = pattern, pb = pb)

  # Identify extant simulations
  if (length(simulations) == 1) simulations <- list.files(simulations, pattern = pattern, full.names = TRUE)
  if (!is.null(missings)) simulations <- simulations[!simulations %in% missings]
  if (!is.null(extincts)) simulations <- simulations[!simulations %in% extincts]

  return (simulations)

}
