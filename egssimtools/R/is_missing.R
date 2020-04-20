#' Is a simulation missing?
#'
#' @param folder Path to the folder
#'
#' @export

is_missing <- function(folder) {

  # A simulation is missing if no data files are available
  length(grep("\\.dat$", list.files(folder))) == 0

}
