#' Is a simulation missing?
#'
#' @param folder Path to the folder
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
#' is_missing(root)
#'
#' }
#'
#' @export

is_missing <- function(folder) {

  # A simulation is missing if no data files are available
  length(grep("\\.dat$", list.files(folder))) == 0

}
