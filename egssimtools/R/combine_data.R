#' Collect data from across simulations
#'
#' Collect specified variables and parameters from across multiple simulation
#' folders
#'
#' @param roots Paths to simulation folders
#' @param type The `type` argument of `read_this`
#' @param id_column Optional name of the simulation identifier column
#' @param archfiles Optional vector of architecture file names, each to be passed with their respective simulation path to the `read_this` function. Must be the same length as `roots` (one architecture for each simulation). The architecture file will be searched for in the simulation folder, so avoid full paths. Only makes sens if `type` is "genome" or "network", ignored otherwise.
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
  roots,
  type = "data",
  id_column = "sim",
  archfiles = NULL,
  ...
) {

  # Pass architecture files along with simulation paths if needed
  if (type %in% c("genome", "network") & !is.null(archfiles)) {

    # Error if not the right number of architecture files provided
    if (length(archfiles) != length(roots)) {
      stop("archfiles must be the same length as roots")
    }

    # Collect data and architectures from all the simulations
    data <- purrr::map2_dfr(
      roots,
      archfiles,
      function(x, y) {
        read_this(x, type, architecture = TRUE, archfile = y, ...)
      },
      .id = id_column
    )

  } else {

    # Read the data and combine
    data <- purrr::map_dfr(roots, read_this, type, ..., .id = id_column)

  }

  return (data)

}
