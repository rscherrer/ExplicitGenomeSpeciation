#' Collect simulation data
#'
#' Collect specified variables and parameters from across multiple simulation folders
#'
#' @param root One of multiple paths to simulation folders or folders into which to recurse to look for simulation folders.
#' @param variables Names of the variables to extract
#' @param by,dupl,parnames,combine,as_numeric Parameters to be passed to `read_data`
#' @param check_extant Whether to check for non-extinct and non-crashed simulation folders
#' @param pattern Optional pattern to look for if simulation folders are searched by recursion
#' @param level Level of recursion. Defaults to 0 for no recursion (hen assumes that `root` is a vector of simulation folder paths).
#' @param verbose Whether to display messages and progress bars
#' @param id_column Optional name of the simulation identifier column
#' @param architecture Whether to read a genetic architecture with the (locus-wise) data
#' @param archfile Name of the architecture file
#' @param parfile Name of the parameter file
#'
#' @return A data frame
#'
#' @examples
#'
#' # Location of the simulation folder
#' root <- "data"
#'
#' collect_sims(
#'   root, c("time", "EI"), pattern = "example", level = 1,
#'   check_extant = FALSE
#' )
#'
#' @export

collect_sims <- function(
  root,
  variables,
  by = 1,
  dupl = 1,
  parnames = NULL,
  combine = FALSE,
  as_numeric = NULL,
  check_extant = TRUE,
  pattern = "sim_",
  level = 0,
  verbose = TRUE,
  id_column = "sim",
  architecture = FALSE,
  archfile = "architecture.txt",
  parfile = "paramlog.txt"
) {

  # Fetch simulation folders
  root <- fetch_dirs(root, pattern = pattern, level = level)

  # Find extant simulations if needed
  if (check_extant) {
    root <- find_extant(root, pattern = pattern, verbose = verbose)
  }

  # Read the data and combine
  if (verbose) message("Reading data...")
  data <- root %>%
    purrr::map_dfr(
      read_data, variables, by, dupl, parnames, combine, as_numeric,
      architecture, archfile, parfile, .id = id_column
    )

  return(data)

}
