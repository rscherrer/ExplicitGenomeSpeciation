#' Collect data from across simulations
#'
#' Collect specified variables and parameters from across multiple simulation
#' folders
#'
#' @param root One of multiple paths to simulation folders or folders into
#' which to recurse to look for simulation folders.
#' @param check_extant Whether to check for non-extinct and non-crashed
#'  simulation folders
#' @param pattern Optional pattern to look for if simulation folders are
#' searched by recursion
#' @param level Level of recursion. Defaults to 0 for no recursion (then assumes
#' that `root` is a vector of simulation folder paths).
#' @param type The `type` argument of `read_this`
#' @param id_column Optional name of the simulation identifier column
#' @param ... Parameters to be passed to `read_this`
#'
#' @return A tibble
#'
#' @examples
#'
#' \dontrun{
#'
#' # Location of the simulation folder
#' root <- "data"
#'
#' # Collect simulation data with some parameters
#' combine_data(
#'   root, variables = c("time", "EI"), pattern = "example", level = 1,
#'   parnames = c("ecosel", "hsymmetry")
#' )
#'
#' combine_data(
#'   root, variables = "genome_Fst", pattern = "example", level = 1,
#'   parnames = c("ecosel", "hsymmetry"), architecture = TRUE, type = "genome"
#' )
#'
#' combine_data(
#'   root, variables = "network_corbreed", pattern = "example", level = 1,
#'   parnames = c("ecosel", "hsymmetry"), architecture = TRUE, type = "network"
#' )
#'
#' }
#'
#' @export

combine_data <- function(
  root,
  check_extant = FALSE,
  pattern = "sim_",
  level = 0,
  type = "data",
  id_column = "sim",
  ...
) {

  # Fetch simulation folders
  if (level > 0) root <- fetch_dirs(root, pattern = pattern, level = level)

  # Find extant simulations if needed
  if (check_extant) root <- find_extant(root)

  # Read the data and combine
  data <- purrr::map_dfr(root, read_this, type, ..., .id = id_column)

  return(data)

}
