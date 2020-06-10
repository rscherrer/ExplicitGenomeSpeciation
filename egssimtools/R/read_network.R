#' Read edge-specific data through time
#'
#' Wrapper around `read_data` to read edge-wise data
#'
#' @param folder Path to the folder
#' @param variables What variables to read (`time` is included by default)
#' @param parnames,combine,as_numeric,parfile Parameters
#' for `read_data`
#' @param architecture Whether to attach genetic architecture data
#' @param archfile Optional name of the genetic architecture file
#' @param nedges Number of edges (automatically guessed if unspecified)
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
#' # Read Fst throughout the genome
#' read_edges(root, "network_corbreed")
#'
#' # Read multiple metrics and attach architecture
#' variables <- paste0("network_", c("corbreed", "corfreq"))
#' read_network(root, variables, architecture = TRUE)
#'
#' }
#'
#' @export

read_network <- function(
  folder,
  variables,
  parnames = NULL,
  combine = FALSE,
  as_numeric = NULL,
  architecture = FALSE,
  archfile = "architecture.txt",
  parfile = "paramlog.txt",
  nedges = NULL
) {

  if (is.null(nedges)) nedges <- guess_nedges(folder, variables[1])
  data <- read_data(
    folder,
    c("time", variables),
    dupl = c(nedges, rep(1, length(variables))),
    parnames = parnames,
    combine = combine,
    as_numeric = as_numeric,
    parfile = parfile
  )
  data$edge <- rep(seq(nedges), nrow(data) / nedges)

  if (architecture) {

    ntimes <- nrow(data) / nedges
    arch <- read_arch_network(folder, archfile, as_list = TRUE)$edges
    arch <- purrr::map_dfr(unique(data$time), ~ dplyr::mutate(arch, time = .x))
    data <- data %>% dplyr::right_join(arch)

  }

  return(data)

}
