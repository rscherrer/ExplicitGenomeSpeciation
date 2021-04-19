#' Is a simulation extinct?
#'
#' Checks for the word "extinct" in the log file
#'
#' @param root Path to the simulation folder
#' @param logfile Optional log file name
#' @param slurm Whether to check a SLURM output file as log file
#'
#' @return A logical
#'
#' @examples
#'
#' \dontrun{
#'
#' is_extinct("data/example_1")
#'
#' }
#'
#' @export

is_extinct <- function(root, logfile = "^log.txt$", slurm = FALSE) {

  if (slurm) logfile <- "^slurm.*out"
  file_id <- grep(logfile, list.files(root))
  if (length(file_id) == 0) stop("log file not found")
  log_file <- list.files(root)[file_id[length(file_id)]]
  log_lines <- readLines(file.path(root, log_file))
  return(length(grep("extinct", log_lines)) > 0)

}

