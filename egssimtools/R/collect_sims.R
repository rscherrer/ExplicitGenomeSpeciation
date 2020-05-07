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
#' @param verbose,pb Whether to display messages and progress bars
#' @param id_column Optional name of the simulation identifier column
#'
#' @return A data frame
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
  pb = TRUE,
  id_column = "sim"
) {

  if (!verbose) pb <- FALSE

  # Fetch simulation folders
  root <- fetch_dirs(root, pattern = pattern, level = level)

  # Find extant simulations if needed
  if (check_extant) {
    root <- find_extant(root, patttern = pattern, verbose = verbose, pb = pb)
  }

  # Read the data and combine
  if (verbose) message("Reading data...")
  data <- root %>%
    purrr::map_dfr(
      read_data, variables, by, dupl, parnames, combine, as_numeric,
      .id = id_column
    )

}
