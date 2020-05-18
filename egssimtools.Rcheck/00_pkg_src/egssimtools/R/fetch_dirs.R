#' Fetch directories within directories
#'
#' @param roots One or several directories to search in
#' @param pattern Optional pattern to look for
#' @param level Recursion level
#'
#' @return A vector of directory full names
#'
#' @export

fetch_dirs <- function(roots, pattern = "^.*$", level = 0) {

  library(tidyverse)

  if (level > 0) for (i in 1:level) roots <- roots %>% list.dirs(recursive = FALSE)
  roots[str_detect(str_replace(roots, "^.*/", ""), pattern)]

}
