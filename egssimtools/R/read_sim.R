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
#' @return A tibble
#'
#' @examples
#'
#' \dontrun{
#'
#' # Location of the simulation folder
#' root <- "data/example_1"
#'
#' # Read the ecological divergence through time
#' read_sim(root, "EI")
#'
#' # Read multiple speciation metrics
#' read_sim(root, c("EI", "RI"))
#'
#' # Read genome-wide metrics
#' read_sim(root, "Fst", by = 3)
#'
#' # Or a combination of multiple metrics with different dimensions
#' read_sim(root, c("EI", "Fst"), by = c(1, 3))
#'
#' # Even locus-specific data in the wide-format
#' read_sim(root, "genome_Fst", by = 300)
#'
#' }
#'
#' @export

read_sim <- function(
  folder,
  variables,
  by = rep(1, length(variables)),
  parnames = NULL,
  combine = FALSE,
  as_numeric = NULL,
  parfile = "paramlog.txt"
) {

  read_data(
    folder,
    variables = c("time", variables),
    by = c(1, by),
    parnames = parnames,
    combine = combine,
    as_numeric = as_numeric,
    parfile = parfile
  )

}
