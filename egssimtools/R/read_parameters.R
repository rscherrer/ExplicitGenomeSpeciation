#' Read parameters
#'
#' @param folder Path to the folder
#' @param parnames Optional vector of parameter names
#' @param filename Optional file name
#' @param combine Whether to paste together compound parameters
#' @param flatten Whether to return a list of single values. If FALSE, some of the elements may be vectors of multiple values (compound parameters).
#'
#' @export

read_parameters <- function(folder, parnames = NULL, filename = "paramlog.txt", combine = FALSE, flatten = FALSE, as_numeric = NULL) {

  library(tidyverse)
  params <- read_paramfile(paste0(folder, '/', filename), parnames)
  if (combine) params <- params %>% map_if(~ length(.x) > 1, ~ paste0(.x, collapse = " "))
  if (flatten) params <- params %>% unlist %>% map(~ .x)
  if (!is.null(as_numeric)) {
    numerics <- as_numeric %>% map(grep, names(params)) %>% unlist
    params[numerics] <- params[numerics] %>% map(as.numeric)
  }
  return (params)

}
