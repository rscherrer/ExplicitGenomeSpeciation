#' Histogram of trait values at a given time
#'
#' Shows the distribution of trait values at a given time.
#'
#' @param folder Path to the folder
#' @param trait What trait to plot? (1) Ecological trait, (2) mating preference and (3) neutral trait.
#' @param t What time point? Defaults to final one
#' @param bins Number of bins in the histogram
#' @param color What color for the histogram?
#' @param is_density Whether to show density
#'
#' @return A ggplot with the basic needed components. You can add extra layers to it (e.g. change the axis labels) outside the function.
#'
#' @export

plot_trait_histogram <- function(folder, trait, t = NULL, bins = 15, color = "lightgrey", is_density = FALSE) {

  library(tidyverse)

  data <- read_individuals(folder, "individual_trait")
  data <- lapply(data, function(data) sapply(data, "[", trait))
  data <- data.frame(time = mrep(as.numeric(names(data)), sapply(data, length)), data = do.call("c", data))

  if (is.null(t)) t <- max(data$time)
  data <- data %>% filter(time == t)

  p <- ggplot(data, aes(x = data))
  if (is_density) p <- p + geom_density(fill = color) else p <- p + geom_histogram(bins = bins, fill = color)
  p <- p +
    theme_bw() +
    xlab("Trait") +
    labs(fill = "Count")
  return (p)

}
