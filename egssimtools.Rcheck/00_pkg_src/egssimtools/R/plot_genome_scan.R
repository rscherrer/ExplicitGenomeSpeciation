#' Genome scan
#'
#' Plot a locus-specific variable across the genome at a given time point
#'
#' @param folder Path to the simulation folder
#' @param variable What locus-specific variable to plot?
#' @param t Time step to plot the scan. Defaults to final time step.
#' @param line Plot lines instead of points
#' @param xaxis What to use for the x-axis: locus index ("locus") or location ("location")?
#' @param colvar Name of a variable to color by. Defaults to NULL for none
#' @param filename Architecture file name. Defaults to architecture.txt
#' @param arch Whether to load a genetic architecture
#'
#' @return A ggplot
#'
#' @export

plot_genome_scan <- function(
  folder,
  variable,
  t = NULL,
  line = FALSE,
  xaxis = "locus",
  colvar = NULL,
  filename = "architecture.txt",
  arch = TRUE,
  col_as_factor = FALSE
) {

  library(tidyverse)

  data <- read_loci(folder, variable)

  if (is.null(t)) t <- last(names(data))
  data <- data[[as.character(t)]]

  data <- data.frame(locus = seq_along(data), y = data)

  # Load genetic architecture if needed
  if (arch) data <- cbind(data, read_genome_architecture(folder, filename))

  # Convert color variable to factor if needed
  if (!is.null(colvar) & col_as_factor) data[, colvar] <- as.factor(data[, colvar])

  p <- ggplot(data, aes(x = get(xaxis), y = y))
  if (line) p <- p + geom_line() else if (is.null(colvar)) p <- p + geom_point() else p <- p + geom_point(aes(color = get(colvar)))
  p <- p +
    theme_bw() +
    ylab(variable) +
    xlab(xaxis)
  if (!is.null(colvar)) p <- p + labs(color = colvar)
  return (p)

}
