#' Collect parameters from multiple simulations
#'
#' @param simulations Paths to the simulations
#' @param parnames optional vector of parameter names (collects all parameters if not specified)
#' @param pb Whether to display a progress bar
#' @param filename Optional file name
#'
#' @return A data frame with parameters in columns and simulations in rows. The parameters are returned as factors.
#'
#' @export

collect_parameters <- function(simulations, parnames = NULL, pb = TRUE, filename = "paramlog.txt") {

  library(pbapply)

  if (pb) thislapply <- pblapply else thislapply <- lapply

  parameters <- thislapply(simulations, read_parameters, parnames, filename)

  # Convert parameter values into factors in a data frame
  parameters <- lapply(parameters, function(parameters) sapply(parameters, function(parameter) paste0(parameter, collapse = " ")))
  parameters <- data.frame(do.call("rbind", parameters))
  rownames(parameters) <- NULL

  return (parameters)

}
