
#' Plot trait values through time
#'
#' Density heatmap of the three trait values across the population.
#'
#' @param data An individual-wise dataset where each trait is a different column.
#' @param bins Number of bins to pass to `ggplot2::geom_bin2d`. Defaults to 30.
#' @param trait_names Vector of trait names as they appear in columns
#'
#' @return A ggplot
#'
#' @export

plot_trait_values <- function(data, bins = 30, trait_names = paste0("individual_trait", 1:3)) {

  data <- tidyr::pivot_longer(data, cols = trait_names, names_to = "trait")

  ggplot2::ggplot(data, ggplot2::aes(x = time, y = value)) +
    ggplot2::geom_bin2d(bins = bins) +
    ggplot2::scale_fill_continuous(type = "viridis") +
    ggplot2::labs(x = "Time (generations)", y = "Trait value", fill = "Count") +
    ggplot2::facet_grid(. ~ trait)

}
