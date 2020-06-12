#' Find extant simulations
#'
#' @param sims A vector of simulation folders
#' @param ... Parameters for `find_extinct`
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
#' roots <- fetch_dirs(root, pattern = "example_", level = 1)
#' find_extant(roots)
#'
#' }
#'
#' @export

find_extant <- function(sims, ...) {

  # Look for missing and extinct simulations
  missings <- find_missing(sims)
  completed <- find_completed(sims)
  extincts <- find_extinct(sims, ...)

  # Identify extant simulations
  if (!is.null(missings)) sims <- sims[!sims %in% missings]
  if (!is.null(completed)) sims <- sims[sims %in% completed]
  if (!is.null(extincts)) sims <- sims[!sims %in% extincts]

  return (sims)

}
