#' Read simulation summary through time
#'
#' Wrapper around `read_data` to read timepoint-wise summary data
#'
#' @param folder Path to the folder
#' @param variables What variables to read (`time` is included by default)
#' @param by Same as the `by` argument of `read_data` (do not count `time` in)
#' @param parnames,combine,as_numeric,parfile Parameters for `read_data`
#'
#' @details The file `time.dat` must be present
#'
#' @return A data frame
#'
#' @examples
#'
#' root <- system.file("extdata", "example_1", package = "egssimtools")
#' read_sim(root, "EI")
#'
#' @export

read_sim <- function(
  folder,
  variables,
  by = 1,
  parnames = NULL,
  combine = FALSE,
  as_numeric = NULL,
  parfile = "paramlog.txt"
) {

  variables <- c("time", variables)
  by <- c(1, by)
  read_data(
    folder, variables, by, parnames = parnames, combine = combine,
    as_numeric = as_numeric, parfile = parfile
  )

}
