#' Fetch directories within directories
#'
#' @param roots One or several directories to search in
#' @param pattern Optional pattern to look for
#' @param level Recursion level
#'
#' @return A vector of directory full names
#'
#' @examples
#'
#' \dontrun{
#'
#' # Location of the simulation folder
#' root <- "data"
#'
#' fetch_dirs(root, pattern = "example", level = 1)
#'
#' }
#'
#' @export

fetch_dirs <- function(roots, pattern = "^.*$", level = 0) {

  if (level > 0) for (i in 1:level) roots <- list.dirs(roots, recursive = FALSE)
  roots[stringr::str_detect(stringr::str_replace(roots, "^.*/", ""), pattern)]

}
