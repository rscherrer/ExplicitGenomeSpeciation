#' Combine parameters from multiple simulations
#'
#' This is essentially a wrapper around `read_parameters`
#'
#' @param root One or multiple paths to simulation folders or folders into
#' which to recurse to look for simulation folders
#' @param parnames,filename,combine,as_numeric Parameters for `read_parameters`
#' @param pattern Optional pattern to look for if simulation folders are
#' searched by recursion
#' @param level Level of recursion. Defaults to 0 for no recursion (then
#' assumes that `root` is a vector of simulation folder paths)
#' @param id_column Optional name of the simulation identifier column
#' @param check_extant Whether to check for non-extinct and non-crashed
#' simulation folders (relevant only if a SLURM output file is present)
#'
#' @return A data frame
#'
#' @examples
#'
#' \dontrun{
#'
#' # Location of the simulation folder
#' root <- "data"
#'
#' collect_param(root, c("ecosel", "hsymmetry"), pattern = "example", level = 1)
#'
#' }
#'
#' @export

collect_param <- function(
  root,
  parnames = NULL,
  filename = "paramlog.txt",
  combine = FALSE,
  as_numeric = NULL,
  pattern = "sim_",
  level = 0,
  id_column = "sim",
  check_extant = FALSE
) {

  # Fetch simulation folders
  if (level > 0) root <- fetch_dirs(root, pattern = pattern, level = level)

  # Find extant simulations if needed
  if (check_extant) root <- find_extant(root)

  # Read the parameters and combine
  purrr::map_dfr(root, function(root) {
    pars <- read_param(
      root, parnames, filename, combine, flatten = TRUE, as_numeric
    )
    do.call("data.frame", pars)
  }, .id = id_column)

}
