#' Did a simulation go extinct?
#'
#' @param folder Path to the folder
#' @param logfile Name of the log file to read
#' @param slurm Whether to look up in the SLURM output file (instead of
#' `logfile`)
#'
#' @details This function looks for the word "extinct" in the log file of the
#' simulation.
#'
#' @return A logical. Returns FALSE and ssues a warning if the log file
#' cannot be found.
#'
#'
#' @examples
#'
#' \dontrun{
#'
#' # Location of the simulation folder
#' root <- "data/example_1"
#'
#' is_extinct(root)
#' is_extinct(root, slurm = TRUE)
#'
#' }
#'
#' @export

is_extinct <- function(folder, logfile = "log.txt", slurm = FALSE) {

  if (slurm) logfile <- "^slurm.*out"

  i <- grep(logfile, list.files(folder))
  if (length(i) == 0) {
    warning(paste("could not find log file for simulation", folder))
    return(FALSE)
  }
  i <- i[length(i)]
  filename <- list.files(folder, full.names = TRUE)[i]
  lines <- read.delim(filename, header = FALSE)[, 1]
  length(grep("extinct", lines)) > 0

}

