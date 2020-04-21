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
  sep = " = "
) {

  library(tidyverse)
  library(pbapply)

  if (!verbose) pb <- FALSE
  if (pb) thislapply <- pblapply else thislapply <- lapply
  if (is.null(colors)) colors <- "black"

  # Look for missing and extinct simulations
  if (verbose) message("Looking for missing simulations...")
  missings <- find_missing(root, pattern = pattern, pb = pb)
  if (verbose) message("Looking for extinct simulations...")
  extincts <- find_extinct(root, pattern = pattern, pb = pb)

  # Identify extant simulations
  simulations <- list.files(root, pattern = pattern, full.names = TRUE)
  if (!is.null(missings)) simulations <- simulations[!simulations %in% missings]
  if (!is.null(extincts)) simulations <- simulations[!simulations %in% extincts]

  # Collect the variable of interest from each simulation
  if (verbose) message("Reading the data...")
  data <- data.frame(do.call("rbind", thislapply(simulations, read_population, variable)))

  # Add time if available
  colnames(data) <- read_time_if(simulations[1])
  timecolumns <- colnames(data)

  # Add simulation identifier
  data$simulation <- factor(rownames(data))

  parnames <- c(facet_rows, facet_cols, color_by)

  # Add parameter values if needed
  if (!is.null(parnames)) {

    if (verbose) message("Reading parameter values...")
    parameters <- thislapply(simulations, read_parameters, parnames)

    # Convert parameter values into factors in a data frame
    parameters <- lapply(parameters, function(parameters) sapply(parameters, function(parameter) paste0(parameter, collapse = " ")))
    parameters <- data.frame(do.call("rbind", parameters))
    rownames(parameters) <- NULL

    # Append to the data
    data <- cbind(data, parameters)

  }

  # Convert from wide to long
  data <- data %>% gather_("time", variable, timecolumns, convert = TRUE)

  ### Choice 1: plot multiple lines ###

  if (!is.null(reverse_order)) {
    if (!all(reverse_order %in% parnames)) stop("Factor to reverse was not provided")
    data[, reverse_order] <- lapply(data.frame(data[, reverse_order]), function(column) factor(column, levels = rev(levels(column))))
  }

  if (verbose) message("Plotting...")

  # Show variable through time
  p <- ggplot(data, aes(x = time, y = get(variable), alpha = simulation))

  # Deal with color
  if (!is.null(color_by)) {

    # Color according to a parameter
    if (color_by_numeric) thisconvert <- function(x) as.numeric(as.character(x)) else thisconvert <- factor
    p <- p + geom_line(aes(color = thisconvert(get(color_by)))) + labs(color = color_by)
    if (color_by_numeric & !is.null(colors)) {
      if (length(colors) != 2) stop("Please provide a lower and upper end for color gradient")
      p <- p + scale_color_gradient(low = colors[1], high = colors[2])
    } else {
      if (length(colors) != nlevels(data[, color_by])) stop("Please provide colors for the different levels of the coloring factor")
      p <- p + scale_color_manual(values = colors)
    }

  } else p <- p + geom_line(color = colors[1])

  p <- p +
    theme_bw() +
    ylab(variable) +
    guides(alpha = FALSE) +
    scale_alpha_manual(values = runif(min = 0.49, max = 0.51, n = nlevels(data$simulation))) # this is a hack to show multiple lines on the same plot

  facets <- c(facet_rows, facet_cols)

  # Use facets if needed
  if (!is.null(facets)) {

    # Setup rows and columns
    if (!is.null(facet_rows)) lhs <- paste(facet_rows, collapse = " + ") else lhs <- "."
    if (!is.null(facet_cols)) rhs <- paste(facet_cols, collapse = " + ") else rhs <- "."

    # Choose what type of facetting to use
    if (facet_wrapped) thisfacet <- facet_wrap else thisfacet <- facet_grid

    # Add custom labels if needed
    if (is.null(facet_prefixes)) facet_prefixes <- facets
    labels <- make_facet_labels(data, facets, facet_prefixes, sep = sep, no_use = !label_facets, simplify = FALSE)

    # Split the plot into facets
    p <- p + thisfacet(formula(paste(lhs, "~", rhs)), labeller = do.call("labeller", labels))

  }

  p

}
