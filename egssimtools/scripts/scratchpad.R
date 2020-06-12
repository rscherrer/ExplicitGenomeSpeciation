rm(list = ls())

library(egssimtools)

roots <- fetch_dirs("data", "example", level = 1)

is_extant <- function(root, slurm = FALSE) {

  if (is_missing(root)) return(FALSE)
  if (slurm & !is_complete(root)) return(FALSE)
  if (is_extinct(root, logfile, slurm)) return(FALSE)
  return(TRUE)

}

find_extant <- function(roots) {

  extant <- purrr::map_lgl(roots, is_extant)
  return(roots[extant])

}

find_extant(roots)
