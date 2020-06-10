#' Did a simulation go extinct?
#'
#' @param folder Path to the folder
#' @param logfile Name of the log file to read
#'
#' @details This function looks for the word "extinct" in the log file of the
#' simulation
#'
#' @return A logical
#'
#' @examples
#'
#' \dontrun{
#'
#' # Location of the simulation folder
#' root <- "data/example_1"
#'
#' is_extinct(root)
#'
#' }
#'
#' @export

is_extinct <- function(folder, logfile = "log.txt") {

  i <- grep(logfile, list.files(folder))
  i <- i[length(i)]
  filename <- list.files(folder, full.names = TRUE)[i]
  lines <- read.delim(filename, header = FALSE)[, 1]
  length(grep("extinct", lines)) > 0

}

