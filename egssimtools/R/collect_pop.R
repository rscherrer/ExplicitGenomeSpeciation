#' Collect individual-wise data from across simulations
#'
#' Collect specified variables and parameters from across multiple simulation
#' folders
#'
#' @param root One of multiple paths to simulation folders or folders into
#' which to recurse to look for simulation folders.
#' @param variables Names of the variables to extract
#' @param by,parnames,combine,as_numeric Parameters to be passed to
#' `read_data`
#' @param check_extant Whether to check for non-extinct and non-crashed
#'  simulation folders
#' @param pattern Optional pattern to look for if simulation folders are
#' searched by recursion
#' @param level Level of recursion. Defaults to 0 for no recursion (then assumes
#' that `root` is a vector of simulation folder paths).
#' @param verbose Whether to display messages and progress bars
#' @param id_column Optional name of the simulation identifier column
#' @param architecture Whether to read a genetic architecture with the
#' (locus-wise) data
#' @param archfile Name of the architecture file
#' @param parfile Name of the parameter file
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
#' collect_pop(
#'   root, "individual_trait", by = 3, pattern = "example", level = 1,
#'   check_extant = FALSE, parnames = c("ecosel", "hsymmetry")
#' )
#'
#' }
#'
#' @export

collect_pop <- function(
  root,
  variables,
  by = rep(1, length(variables)),
  parnames = NULL,
  combine = FALSE,
  as_numeric = NULL,
  check_extant = TRUE,
  pattern = "sim_",
  level = 0,
  verbose = TRUE,
  id_column = "sim",
  parfile = "paramlog.txt"
) {

  collect_data(
    root,
    variables = c("time", variables),
    by = c(1, by),
    dupl = c("population_size", as.list(rep(1, length(variables)))),
    parnames = parnames,
    combine = combine,
    as_numeric = as_numeric,
    check_extant = check_extant,
    pattern = pattern,
    level = level,
    verbose = verbose,
    id_column = id_column,
    parfile = parfile
  )

}
