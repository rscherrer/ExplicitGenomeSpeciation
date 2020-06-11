#' Find extant simulations
#'
#' @param sims Either path to a root directory containing simulation folders,
#' or a vector of simulation folders
#' @param pattern Pattern defining the simulation folders to look into
#' @param verbose Whether to display messages
#' @param level Recursion level
#' @param slurm Parameter for `find_extinct`
#'
#' @return A vector of names of the simulation that did not go extinct
#'
#' @examples
#'
#' \dontrun{
#'
#' # Location of the simulation folder
#' root <- "data"
#'
#' find_extant(root, pattern = "example_", level = 1)
#'
#' }
#'
#' @export

find_extant <- function(
  sims,
  pattern = "^sim_",
  verbose = TRUE,
  level = 0,
  slurm = FALSE
) {

  if (verbose) message("Looking for extant simulations...")
  if (level > 0) sims <- fetch_dirs(
    sims, pattern = pattern, level = level
  )

  # Look for missing and extinct simulations
  missings <- find_missing(sims, verbose = verbose)
  completed <- find_completed(sims, verbose = verbose)
  extincts <- find_extinct(sims, verbose = verbose, slurm = slurm)

  # Identify extant simulations
  if (!is.null(missings)) sims <- sims[!sims %in% missings]
  if (!is.null(completed)) sims <- sims[sims %in% completed]
  if (!is.null(extincts)) sims <- sims[!sims %in% extincts]

  return (sims)

}
