#' Collect parameters from multiple simulations
#'
#' @param simulations Either path to a root directory containing simulation folders, or a vector of simulation folders
#' @param parnames Optional vector of parameter names (collects all parameters if not specified)
#' @param pattern Optional pattern characteristic of simulation folders. Defaults to starting with "sim_".
#' @param verbose Whether to display messages
#' @param pb Whether to display a progress bar
#' @param filename Optional file name
#' @param check_extant Whether to filter extant simulations. Defaults to TRUE. Set to FALSE if e.g. you supplied a vector of simulations you know did not go missing or extinct.
#' @param to_numeric Which parameters to convert into numeric. Because all parameters may not be numbers, they are read as factor by default.
#' @param add_sim Whether to add simulation identifiers
#' @param as_address Whether to use the path to the simulation folders as simulation identifiers in the resulting data frame. Defaults to FALSE, but useful to retrieve a particular simulation location
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
  check_extant = TRUE,
  to_numeric = NULL,
  add_sim = FALSE,
  as_address = FALSE
) {

  library(pbapply)
  library(tidyverse)

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

  if (add_sim) parameters$simulation <- factor(seq_along(simulations))
  if (add_sim & as_address) parameters$simulation <- factor(simulations)

  # Convert specific parameters into numeric
  if (!is.null(to_numeric)) {
    parnames <- colnames(parameters)
    if (!all(to_numeric %in% colnames(parameters))) stop("Unknown parameter to convert into numeric")
    parameters <- parameters %>% mutate_at(to_numeric, function(x) as.numeric(as.character(x)))
  }

  return (parameters)

}
