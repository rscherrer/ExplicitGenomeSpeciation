
#' Plot trait values through time
#'
#' Density heatmap of the three trait values across the population.
#'
#' @param data An individual-wise dataset where each trait is a different column.
#'
#' @return A ggplot
#'
#' @export

plot_trait_values <- function(data) {

  data <- tidyr::pivot_longer(data, cols = paste0("individual_trait", 1:3), names_to = "trait")
  data <- data %>% dplyr::mutate(trait = stringr::str_remove(trait, "individual_trait"))

  ggsim::facettize(

    ggplot2::ggplot(data, ggplot2::aes(x = time, y = value)) +
      ggplot2::geom_bin2d() +
      ggplot2::scale_fill_continuous(type = "viridis") +
      ggplot2::labs(x = "Time (generations)", y = "Trait value", fill = "Count"),

    cols = "trait",
    prepend = "Trait "

  )

}
