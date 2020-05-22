#' Plot densities from genome data
#'
#' Can plot the density of a variable across the genome in various ways
#'
#' @param folder The simulation folder
#' @param variable What variable to plot?
#' @param t What time point to plot? Defaults to the last one
#' @param filename Name of the architecture file, if needed
#' @param arch Whether to read a genetic architecture
#' @param graph List of named graphical parameters to be passed to `ggdensityplot`
#'
#' @return A ggplot
#'
#' @export

plot_genome_distribution <- function(
  folder,
  variable,
  t = NULL,
  filename = "architecture.txt",
  arch = TRUE,
  graph = NULL
) {

  library(tidyverse)

  data <- read_loci(folder, variable)

  if (is.null(t)) t <- last(names(data))
  data <- data[[as.character(t)]]

  data <- data.frame(locus = seq_along(data), y = data)
  colnames(data)[colnames(data) == "y"] <- variable

  # Load genetic architecture if needed
  if (arch) data <- cbind(data, read_genome_architecture(folder, filename))

  do.call("ggdensityplot", c(list(data), variable, graph))

}
