#' Read locus-specific data
#'
#' @param folder Path to the folder
#' @param variable What to read
#'
#' @export

read_loci <- function(folder, variable) {

  t <- read_time(folder)
  data <- read_data(folder, variable)
  nloci <- length(data) / length(t)
  timepoints <- rep(t, each = nloci)
  split(data, timepoints)

}
