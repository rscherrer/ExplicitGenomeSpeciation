#' Read locus-specific data through time
#'
#' Wrapper around `read_data` to read locus-wise data
#'
#' @param folder Path to the folder
#' @param variables What variables to read (`time` is included by default)
#' @param parnames,combine,as_numeric,architecture,archfile,parfile Parameters for `read_data`
#' @param nloci Number of loci (automatically guessed if unspecified)
#'
#' @details The file `time.dat` must be present
#'
#' @return A data frame
#'
#' @examples
#'
#' # Location of the simulation folder
#' root <- "egsimtools/data/example_1"
#'
#' # Read Fst throughout the genome
#' read_genome(root, "genome_Fst")
#'
#' # Read multiple metrics and attach architecture
#' read_genome(root, c("genome_Fst", "genome_Cst"), architecture = TRUE)
#'
#' @export

read_genome <- function(
  folder,
  variables,
  parnames = NULL,
  combine = FALSE,
  as_numeric = NULL,
  architecture = FALSE,
  archfile = "architecture.txt",
  parfile = "paramlog.txt",
  nloci = NULL
) {

  if (is.null(nloci)) nloci <- guess_nloci(folder, variables[1])

  data <- read_data(
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

  if (!"locus" %in% colnames(data)) {
    data$locus <- rep(seq(nloci), nrow(data) / nloci)
  }

  return(data)

}
