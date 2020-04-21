#' Collect data across simulations
#'
#' Read in population-wide data and parameters from multiple extant simulations
#'
#' @param root Where the simulation folders are
#' @param variable What variable to plot
#' @param pattern Optional pattern characteristic of simulation folders. Defaults to starting with "sim_".
#' @param pb Whether to display progress bars
#' @param verbose Whether to display messages
#' @param reverse_order Optional name of a parameter for which to reverse the order of the facets (e.g. to put highest values on the top rows)
#' @param add_summaries Optional function to add extra columns to the data prior to plotting. This can be useful for mapping more complex values than those read from parameter files, to aesthetics such as facets or colors (e.g. color the lines by value of the mean of the variable over the final few timepoint). The function will be called with a single argument, the long data frame to which ggplot is applied, and must return a table with columns to append to that long data frame. Those added columns can be specified and used in facet_rows, facet_cols and color_by, just as any parameter. Make sure that potential extra arguments are passed by default or within the function body (e.g. the time points over which to measure the mean of the variable). Its output will be appended to the long data frame using cbind(), so make sure that it returns a table with new, summary variables as columns, and the right number of rows. The long table taken as input has at least the columns "simulation" (factor), "time" (integer), the variable to plot and any optional parameters read from parameter files that are requested in facet_rows, facet_cols or color_by. As an example, "add_summaries = function(data) data %>% group_by(simulation) %>% mutate(x = last(RI)) %>% ungroup() %>% select(x)" will add a column named "x" containing the last value of variable RI for each simulation, and assumes that RI is the variable to be plotted here. Note that the "plot_simulations" function loads the tidyverse, so no need to load any of it in the function passed to this argument.
#'
#' @return A data frame in the long format, ready for plotting. It contains columns with the variable of interest, time, simulation identifier, parameters if requested and optional extra variables (see add_summaries).
#'
#' @export


collect_simulations <- function(
  root,
  variable,
  parnames = NULL,
  pattern = "^sim_",
  verbose = TRUE,
  pb = TRUE,
  add_summaries = NULL,
  reverse_order = NULL
) {

  library(pbapply)
  library(tidyverse)

  if (!verbose) pb <- FALSE
  if (pb) thislapply <- pblapply else thislapply <- lapply

  # Find the simulations that completed
  simulations <- find_extant(root, pattern = pattern, verbose = verbose, pb = pb)

  # Collect the variable of interest from each simulation
  if (verbose) message("Reading the data...")
  data <- data.frame(do.call("rbind", thislapply(simulations, read_population, variable)))

  # Add time if available
  colnames(data) <- read_time_if(simulations[1])
  timecolumns <- colnames(data)

  # Add simulation identifier
  data$simulation <- factor(rownames(data))

  # Add parameter values if needed
  if (!is.null(parnames)) {
    if (verbose) message("Reading parameter values...")
    data <- cbind(data, collect_parameters(simulations, parnames))
  }

  # Convert from wide to long
  data <- data %>% gather_("time", variable, timecolumns, convert = TRUE)

  # Add summaries if needed
  if (!is.null(add_summaries)) data <- cbind(data, add_summaries(data))

  # Reverse the order of some factors if needed (e.g. nicer facet or legend plotting)
  if (!is.null(reverse_order)) {
    if (!all(reverse_order %in% parnames)) stop("Factor to reverse was not provided")
    data[, reverse_order] <- lapply(data.frame(data[, reverse_order]), function(column) factor(column, levels = rev(levels(column))))
  }

  return (data)

}
