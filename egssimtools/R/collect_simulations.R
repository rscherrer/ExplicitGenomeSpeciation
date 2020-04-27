#' Collect data across simulations
#'
#' Read in population-wide data and parameters from multiple extant simulations
#'
#' @param simulations Either path to a root directory containing simulation folders, or a vector of simulation folders
#' @param variables What variables to collect
#' @param parnames Optional vector of parameter names (collects all parameters if not specified)
#' @param to_numeric Which parameters to convert into numeric. Because all parameters may not be numbers, they are read as factor by default.
#' @param pattern Optional pattern characteristic of simulation folders. Defaults to starting with "sim_".
#' @param pb Whether to display progress bars
#' @param verbose Whether to display messages
#' @param reverse_order Optional name of a parameter for which to reverse the order of the levels (e.g. to put highest values on the top rows)
#' @param add_summaries Optional function to add extra columns to the data prior to plotting. This can be useful for mapping more complex values than those read from parameter files, to aesthetics such as facets or colors (e.g. color the lines by value of the mean of the variable over the final few timepoint). The function will be called with a single argument, the long data frame to which ggplot is applied, and must return a table with columns to append to that long data frame. Those added columns can be specified and used in facet_rows, facet_cols and color_by, just as any parameter. Make sure that potential extra arguments are passed by default or within the function body (e.g. the time points over which to measure the mean of the variable). Its output will be appended to the long data frame using cbind(), so make sure that it returns a table with new, summary variables as columns, and the right number of rows. The long table taken as input has at least the columns "simulation" (factor), "time" (integer), the variable to plot and any optional parameters read from parameter files that are requested in facet_rows, facet_cols or color_by. As an example, "add_summaries = function(data) data %>% group_by(simulation) %>% mutate(x = last(RI)) %>% ungroup() %>% select(x)" will add a column named "x" containing the last value of variable RI for each simulation, and assumes that RI is the variable to be plotted here. Note that the "plot_simulations" function loads the tidyverse, so no need to load any of it in the function passed to this argument.
#' @param filters Optional strings to be parsed into a call to the dplyr::filter function, allowing various rules to filter the data. For example, use filters = c("ecosel == 1", "hsymmetry %in% c(0, 1)") to only keep simulations where ecosel is 1 and hsymmetry is either 0 or 1. Those parsed expressions must evaluate to logicals when the function is called.
#' @param check_extant Whether to filter extant simulations. Defaults to TRUE. Set to FALSE if e.g. you supplied a vector of simulations you know did not go missing or extinct.
#' @param as_address Whether to use the path to the simulation folders as simulation identifiers in the resulting data frame. Defaults to FALSE, but useful to retrieve a particular simulation location
#'
#' @details At least root or simulations must be provided.
#'
#' @return A data frame in the long format. It contains columns with the variables of interest, time, simulation identifier, parameters if requested and optional extra variables (see add_summaries).
#'
#' @export


collect_simulations <- function(
  simulations,
  variables,
  parnames = NULL,
  to_numeric = NULL,
  pattern = "^sim_",
  verbose = TRUE,
  pb = TRUE,
  add_summaries = NULL,
  reverse_order = NULL,
  filters = NULL,
  check_extant = TRUE,
  as_address = FALSE
) {

  library(pbapply)
  library(tidyverse)

  if (!verbose) pb <- FALSE
  if (pb) thislapply <- pblapply else thislapply <- lapply

  # Find the simulations that completed
  if (length(simulations) == 1) simulations <- list.files(simulations, pattern = pattern, full.names = TRUE)
  if (check_extant) simulations <- find_extant(simulations, verbose = verbose, pb = pb)

  # Collect the variable of interest from each simulation
  if (verbose) message("Reading the data...")
  data <- data.frame(do.call("rbind", thislapply(simulations, read_population, variables[1])))

  # Add time if available
  colnames(data) <- read_time_if(simulations[1])
  timecolumns <- colnames(data)

  # Add simulation identifier
  data$simulation <- factor(rownames(data))
  if (as_address) data$simulation <- factor(simulations)

  # Add parameter values if needed
  if (!is.null(parnames)) data <- cbind(collect_parameters(simulations, parnames, pattern = pattern, verbose = verbose, check_extant = FALSE, to_numeric = to_numeric), data)

  # Convert from wide to long
  data <- data %>% gather_("time", variables[1], timecolumns, convert = TRUE)

  # Add optional extra variables
  if (length(variables) > 1) {

    if (verbose) message("Reading extra data...")

    xdata <- do.call("cbind", lapply(variables[2:length(variables)], function(variable) {

      xdata <- data.frame(do.call("rbind", thislapply(simulations, read_population, variable)))
      colnames(xdata) <- timecolumns
      xdata <- xdata %>% gather_("time", variable, timecolumns) %>% select(all_of(variable))

    }))

    data <- cbind(data, xdata)
  }

  # Add summaries if needed
  if (!is.null(add_summaries)) data <- cbind(data, add_summaries(data))

  # Reverse the order of some factors if needed (e.g. nicer facet or legend plotting)
  if (!is.null(reverse_order)) {
    if (!all(reverse_order %in% parnames)) stop("Factor to reverse was not provided")
    data[, reverse_order] <- lapply(data.frame(data[, reverse_order]), function(column) factor(column, levels = rev(levels(column))))
  }

  # Filter if needed
  if (!is.null(filters)) {
    eval(parse(text = sprintf("thisfilter <- function(data) { filter(data, %s) }", paste0(filters, collapse = ", "))))
    data <- thisfilter(data)
  }

  return (data)

}

