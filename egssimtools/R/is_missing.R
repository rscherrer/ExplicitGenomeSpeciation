#' Is a simulation missing?
#'
#' Checks for the presence of `.dat` files
#'
#' @param root Path to the simulation folder
#'
#' @return A logical
#'
#' @examples
#'
#' \dontrun{
#'
#' is_missing("data/example_1")
#'
#' }
#'
#' @export

is_missing <- function(root) {

  if (length(grep(".dat$", list.files(root))) == 0) return(TRUE)
  return(FALSE)

}
