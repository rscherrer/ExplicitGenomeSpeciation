#' Read the SLURM status of a simulation
#'
#' Checks the slurm output file
#'
#' @param root Path to the simulation folder
#'
#' @return A character
#'
#' @examples
#'
#' \dontrun{
#'
#' read_slurm_status("data/example_1") # should give an error because examples
#' # were not run using SLURM
#'
#' }
#'
#' @export

read_slurm_status <- function(root) {

  file_id <- grep("^slurm.*out$", list.files(root))
  if (length(file_id) == 0) stop("SLURM output file not found")
  slurm_file <- list.files(root)[file_id[length(file_id)]]
  slurm_lines <- readLines(file.path(root, slurm_file))
  i <- grep("^State", slurm_lines)
  if (length(i) == 0) {
    stop(paste("row starting with 'State' not found in", slurm_file))
  }
  status_line <- slurm_lines[i]
  status <- gsub("^.*\\: ", "", status_line)
  return(status)

}
