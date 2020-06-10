#' Read parameters
#'
#' @param folder Path to the folder
#' @param parnames Optional vector of parameter names
#' @param filename Optional file name
#' @param combine Whether to paste together compound parameters
#' @param flatten Whether to return a list of single values. If FALSE, some of
#' the elements may be vectors of multiple values (compound parameters).
#' @param as_numeric Which parameters to convert to numeric?
#'
#' @return A list of parameter values
#'
#' @examples
#'
#' \dontrun{
#'
#' # Location of the simulation folder
#' root <- "data/example_1"
#'
#' read_param(root)
#'
#' }
#'
#' @export

read_param <- function(
  folder,
  parnames = NULL,
  filename = "paramlog.txt",
  combine = FALSE,
  flatten = FALSE,
  as_numeric = NULL
) {

  params <- read_param_file(paste0(folder, '/', filename), parnames)
  if (combine) params <- purrr::map_if(
    params, ~ length(.x) > 1, ~ paste0(.x, collapse = " ")
  )
  if (flatten) params <- params %>% unlist %>% purrr::map(~ .x)
  if (!is.null(as_numeric)) {
    numerics <- as_numeric %>% purrr::map(grep, names(params)) %>% unlist
    params[numerics] <- params[numerics] %>% purrr::map(as.numeric)
  }
  return(params)

}
