#' Collect parameters from multiple simulations
#'
#' @param simulations Either path to a root directory containing simulation folders, or a vector of simulation folders
#' @param parnames optional vector of parameter names (collects all parameters if not specified)
#' @param pattern Optional pattern characteristic of simulation folders. Defaults to starting with "sim_".
#' @param verbose Whether to display messages
#' @param pb Whether to display a progress bar
#' @param filename Optional file name
#' @param check_extant Whether to filter extant simulations. Defaults to TRUE. Set to FALSE if e.g. you supplied a vector of simulations you know did not go missing or extinct.
#'
#' @return A data frame with parameters in columns and simulations in rows. The parameters are returned as factors.
#'
#' @export

collect_parameters <- function(
  simulations,
  parnames = NULL,
  pattern = "^sim_",
  verbose = TRUE,
  pb = TRUE,
  filename = "paramlog.txt",
  check_extant = TRUE
) {

  library(pbapply)

  if (!verbose) pb <- FALSE
  if (pb) thislapply <- pblapply else thislapply <- lapply

  if (verbose) message("Reading parameter values...")
  if (length(simulations) == 1) simulations <- list.files(simulations, pattern = pattern, full.names = TRUE)
  if (check_extant) simulations <- find_extant(simulations, verbose = verbose, pb = pb)
  parameters <- thislapply(simulations, read_parameters, parnames, filename)

  # Convert parameter values into factors in a data frame
  parameters <- lapply(parameters, function(parameters) sapply(parameters, function(parameter) paste0(parameter, collapse = " ")))
  parameters <- data.frame(do.call("rbind", parameters))
  rownames(parameters) <- NULL

  return (parameters)

}
