#' Heatmap of trait values through time
#'
#' @param filename Path to the folder
#' @param trait What trait to plot? (1) Ecological trait, (2) mating preference and (3) neutral trait.
#' @param bins Number of bins in the heatmap
#'
#' @return A ggplot with the basic needed components. You can add extra layers to it (e.g. change the axis labels) outside the function.
#'
#' @export

plot_trait_heatmap <- function(folder, trait, bins = 15) {

  library(tidyverse)

  # Plot a simulation through time
  data <- read_individuals(folder, "individual_trait")
  data <- lapply(data, function(data) sapply(data, "[", trait))
  data <- data.frame(time = mrep(as.numeric(names(data)), sapply(data, length)), data = do.call("c", data))

  ggplot(data, aes(x = time, y = data)) +
    geom_bin2d(bins = bins) +
    theme_bw() +
    scale_fill_continuous(type = "viridis") +
    xlab("Time (generations)") +
    labs(fill = "Count")

}
