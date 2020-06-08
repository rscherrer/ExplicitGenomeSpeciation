#' Is a simulation missing?
#'
#' @param folder Path to the folder
#'
#' @return A logical
#'
#' @examples
#'
#' root <- system.file("extdata", "example_1", package = "egssimtools")
#' is_missing(root)
#'
#' @export

is_missing <- function(folder) {

  # A simulation is missing if no data files are available
  length(grep("\\.dat$", list.files(folder))) == 0

}
