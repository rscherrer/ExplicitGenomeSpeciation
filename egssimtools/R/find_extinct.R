#' Find extinct simulations
#'
#' @param sims Either path to a root directory containing simulation folders,
#' or a vector of simulation folders
#' @param pattern Pattern defining the simulation folders to look into
#' @param verbose Whether to display messages
#' @param level Recursion level
#'
#' @return A character vector of extinct simulations, or NULL if no simulation
#'   is extinct
#'
#' @examples
#'
#' \dontrun{
#'
#' # Location of the simulation folder
#' root <- "data"
#'
#' find_extinct(root, pattern = "example_", level = 1)
#'
#' }
#'
#' @export

# How many simulations did go extinct?
find_extinct <- function(
  sims,
  pattern = "^sim_",
  verbose = TRUE,
  level = 0
) {

  if (verbose) message("Looking for extinct simulations...")
  if (level > 0) sims <- fetch_dirs(
    sims, pattern = pattern, level = level
  )

  # Check each one for extinction
  extincts <- sapply(sims, is_extinct)

  # Return the names of the extinct folders
  if (any(extincts)) return (sims[extincts])
  return (NULL)

}

