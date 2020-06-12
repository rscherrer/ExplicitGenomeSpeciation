#' Collect SLURM exit status
#'
#' @param sims A vector of simulation folders
#'
#' @return A character vector
#'
#' @examples
#'
#' \dontrun{
#'
#' # Location of the simulation folder
#' root <- "data"
#'
#' roots <- fetch_dirs(root, pattern = "example", level = 1)
#'
#' # Should return not found with the example data because these simulations
#' # were run locally, not on SLURM
#' collect_status(roots)
#'
#' }
#'
#' @export

collect_status <- function(sims) {

  return (purrr::map_chr(sims, get_status))

}
