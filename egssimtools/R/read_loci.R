#' Read locus-specific data through time
#'
#' Wrapper around `read_data` to read locus-wise data
#'
#' @param folder Path to the folder
#' @param variables What variables to read (`time` is included by default)
#' @param parnames,combine,as_numeric,architecture,archfile,parfile Parameters for `read_data`
#'
#' @details The file `time.dat` must be present
#'
#' @return A data frame
#'
#' @examples
#'
#' root <- system.file("extdata", "example_1", package = "egssimtools")
#' read_loci(root, "genome_Fst")
#'
#' @export

read_loci <- function(
  folder,
  variables,
  parnames = NULL,
  combine = FALSE,
  as_numeric = NULL,
  architecture = FALSE,
  archfile = "architecture.txt",
  parfile = "paramlog.txt"
) {

  nloci <- guess_nloci(folder, variables[1])

  read_data(
    folder,
    c("time", variables),
    by = rep(1, length(variables) + 1),
    dupl = c(nloci, rep(1, length(variables))),
    parnames = parnames,
    combine = combine,
    as_numeric = as_numeric,
    architecture = architecture,
    archfile = archfile,
    parfile = parfile
  )

}
