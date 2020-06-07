#' Read individual data through time
#'
#' Wrapper around `read_data` to read individual-wise data
#'
#' @param folder Path to the folder
#' @param variables What variables to read (`time` is included by default)
#' @param by Same as the `by` argument of `read_data` (do not count `time` in)
#' @param parnames,combine,as_numeric,parfile Parameters for `read_data`
#'
#' @details The files `time.dat` and `population_size.dat` must be present
#'
#' @return A data frame
#'
#' @export

read_indiv <- function(
  folder,
  variables,
  by = 3,
  parnames = NULL,
  combine = FALSE,
  as_numeric = NULL,
  parfile = "paramlog.txt"
) {

  read_data(
    folder,
    c("time", variables),
    by = c(1, by),
    dupl = list("population_size", rep(1, length(variables))),
    parnames = parnames,
    combine = combine,
    as_numeric = as_numeric,
    parfile = parfile
  )

}
