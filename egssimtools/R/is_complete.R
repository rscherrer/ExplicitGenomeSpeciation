#' Is a simulation completed?
#'
#' Checks for the word "COMPLETED" in the SLURM output
#'
#' @param root Path to the simulation folder
#'
#' @return A logical
#'
#' @examples
#'
#' \dontrun{
#'
#' # Should return an error because examples were not run using SLURM
#' is_complete("data/example_1")
#'
#' }
#'
#' @export

is_complete <- function(root) {

  if (length(grep("COMPLETED", read_slurm_status(root))) > 0) return(TRUE)
  return(FALSE)

}
