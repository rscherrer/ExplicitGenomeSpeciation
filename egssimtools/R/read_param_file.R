#' Read a parameter file
#'
#' @param filename Path to the parameter file
#' @param parnames Optional vector of parameter names
#'
#' @return A list of parameter values
#'
#' @examples
#'
#' \dontrun{
#'
#' # Location of the simulation folder
#' root <- "data/example_1"
#'
#' f <- file.path(root, "parameters.txt")
#' read_param_file(f)
#'
#' }
#'
#' @export

read_param_file <- function(filename, parnames = NULL) {

  # Give it the path to the parameter file
  # And optionally the parameters you want it to read

  paramlist <- readLines(filename)

  # Extract values of the parameters
  parameters <- lapply(paramlist, function(p) {

    p <- strsplit(as.character(p), " ")[[1]]

    # Get parameter only if required
    if (!is.null(parnames) & !p[1] %in% parnames) return (NULL)
    return (list(p[1], p[-1]))

  })

  parameters <- parameters[sapply(parameters, function(x) !is.null(x))]
  names(parameters) <- sapply(parameters, "[[", 1)
  parameters <- lapply(parameters, "[[", 2)

  return (parameters)

}
