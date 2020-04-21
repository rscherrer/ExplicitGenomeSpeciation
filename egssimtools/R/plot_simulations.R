#' Plot multiple simulations through time
#'
#' This function allows to easily loop through many simulation folders, extract data from them and plot them synthetically. Some population-wide variables are plotted through time for each simulation.
#'
#' @param root Where the simulation folders are
#' @param variable What variable to plot
#' @param pattern Optional pattern characteristic of simulation folders. Defaults to starting with "sim_".
#' @param color_by Optional parameter to use to color the lines
#' @param color_by_numeric Whether the color parameter must be treated as numeric and be mapped to a color gradient. Defaults to TRUE. If FALSE, it is treated as a factor.
#' @param colors Optional colors of the lines. If color_by is NULL, the color that all the lines must have. If color_by is defined and color_by_numeric is TRUE, the lower and the upper end of the color gradient. If color_by_numeric is FALSE, a set of colors for each level of the color factor. If not defined, default ggplot color(s) will be used.
#' @param pb Whether to display progress bars
#' @param verbose Whether to display messages
#' @param facet_rows Optional parameter name(s) to facet the plot by rows
#' @param facet_cols Optional parameter name(s) to facet the plot by columns
#' @param facet_wrapped Whether to automatically fill the facets by row with one-dimensional array of parameter combinations. If TRUE, combines the parameters in facet_rows and facet_cols into a single array. Ignored if none of those is specified.
#' @param reverse_order Optional name of a parameter for which to reverse the order of the facets (e.g. to put highest values on the top rows)
#' @param label_facets Whether to add custom labels to the facets
#' @param facet_prefixes Optional prefixes to add to the facet labels of each parameter. Must have one value per facetting parameter. Ignored if label_facets is FALSE. If not specified and label_facets is TRUE, the names of the parameters are used as prefixes.
#' @param sep Optional separator to use when adding prefixes. Defaults to " = ".
#' @param add_summaries Optional function to add extra columns to the data prior to plotting. This can be useful for mapping more complex values than those read from parameter files, to aesthetics such as facets or colors (e.g. color the lines by value of the mean of the variable over the final few timepoint). The function will be called with a single argument, the long data frame to which ggplot is applied, and must return a table with columns to append to that long data frame. Those added columns can be specified and used in facet_rows, facet_cols and color_by, just as any parameter. Make sure that potential extra arguments are passed by default or within the function body (e.g. the time points over which to measure the mean of the variable). Its output will be appended to the long data frame using cbind(), so make sure that it returns a table with new, summary variables as columns, and the right number of rows. The long table taken as input has at least the columns "simulation" (factor), "time" (integer), the variable to plot and any optional parameters read from parameter files that are requested in facet_rows, facet_cols or color_by. As an example, "add_summaries = function(data) data %>% group_by(simulation) %>% mutate(x = last(RI)) %>% ungroup() %>% select(x)" will add a column named "x" containing the last value of variable RI for each simulation, and assumes that RI is the variable to be plotted here. Note that the "plot_simulations" function loads the tidyverse, so no need to load any of it in the function passed to this argument.
#' @param filters Optional strings to be parsed into a call to the dplyr::filter function, allowing various rules to filter the data. For example, use filters = c("ecosel == 1", "hsymmetry %in% c(0, 1)") to only keep simulations where ecosel is 1 and hsymmetry is either 0 or 1. Those parsed expressions must evaluate to logicals when the function is called.
#'
#' @return A plot showing a variable through time for multiple simulations.
#'
#' @note The lines have different transparency levels, all very close to 0.5. This is a hack to make sure that all simulations can efficiently be plotted together on the same plot.
#'
#' @export

plot_simulations <- function(
  root,
  variable,
  pattern = "^sim_",
  color_by = NULL,
  color_by_numeric = TRUE,
  colors = NULL,
  pb = TRUE,
  verbose = TRUE,
  facet_rows = NULL,
  facet_cols = NULL,
  facet_wrapped = FALSE,
  reverse_order = NULL,
  label_facets = FALSE,
  facet_prefixes = NULL,
  sep = " = ",
  add_summaries = NULL,
  filters = NULL
) {

  library(tidyverse)
  library(pbapply)

  parnames <- c(facet_rows, facet_cols, color_by)

  # Collect simulation data in the long format
  data <- collect_simulations(root, variable, parnames, pattern, verbose, pb, add_summaries, reverse_order, filters)

  if (verbose) message("Plotting...")
  p <- ggplot(data, aes(x = time, y = get(variable), alpha = simulation))

  # Deal with color
  p <- add_geom_line_colored(p, color_by, color_by_numeric, colors)

  # Add general plot components
  p <- p +
    theme_bw() +
    ylab(variable) +
    guides(alpha = FALSE) +
    scale_alpha_manual(values = runif(min = 0.49, max = 0.51, n = nlevels(data$simulation)))

  # Add facets if needed (the function handle no facet too)
  facettize(p, facet_rows, facet_cols, facet_wrapped, label_facets, facet_prefixes, sep)

}
