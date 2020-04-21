#' 2D heatmap of trait densities
#'
#' @param filename Path to the folder
#' @param traits Vector of two traits to plot: (1) ecological trait, (2) mating preference and (3) neutral trait.
#' @param t What time point to show. All time points together if NULL.
#' @param bins Number of bins in the heatmap
#'
#' @return A ggplot with the basic needed components. You can add extra layers to it (e.g. change the axis labels) outside the function.
#'
#' @export

plot_trait_heatmap2D <- function(folder, traits = c(1, 2), t = NULL, bins = 15) {

  library(tidyverse)

  # Plot a simulation through time
  data <- read_individuals(folder, "individual_trait")
  data <- lapply(traits, function(trait) {
    lapply(data, function(data) sapply(data, "[", trait))
  })
  data <- data.frame(
    time = mrep(as.numeric(names(data[[1]])), sapply(data[[1]], length)),
    do.call("cbind", lapply(data, function(data) do.call("c", data)))
  )

  if (!is.null(t)) data <- data %>% filter(time == t)

  ggplot(data, aes(x = X1, y = X2)) +
    geom_bin2d(bins = bins) +
    theme_bw() +
    scale_fill_continuous(type = "viridis") +
    labs(fill = "Count")

}
