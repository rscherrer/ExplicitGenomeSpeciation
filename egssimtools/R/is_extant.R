#' Is a simulation extant?
#'
#' Checks that the simulation is non-missing, succesfully completed and did
#' not go extinct.
#'
#' @param root Path to the simulation folder
#' @param logfile Optional log file name
#' @param slurm Whether to check SLURM output files
#'
#' @return A logical
#'
#' @examples
#'
#' \dontrun{
#'
#' is_extant("data/example_1", slurm = FALSE)
#'
#' }
#'
#' @export

is_extant <- function(root, logfile = "^log.txt$", slurm = FALSE) {

  if (is_missing(root)) return(FALSE)
  if (slurm) if(!is_complete(root)) return(FALSE)
  return(!is_extinct(root, logfile = logfile, slurm = slurm))

}

