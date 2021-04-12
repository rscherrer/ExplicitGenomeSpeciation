#' Add chromosomes to plot
#'
#' Add chromosomes as lines to a genome scan.
#'
#' @param p A ggplot showing a genome scan (see `?plot_genome_scan`).
#' @param y Vertical location of the chromosomes on the plot.
#' @param offset Small number controlling the distance between the chromosomes.
#' @param size Width of the bars.
#' @param color Color of the bars.
#' @param unit Whether the horizontal axis is "locus" or "location" (see `?plot_genome_scan`).
#'
#' @return A ggplot
#'
#' @export

add_chromosomes <- function(
  p, y = -0.05, offset = 0.01, size = 2, color = "black", unit = "location"
) {

  # Extract the limits of the chromosomes
  chromosomes <- p$data %>%
    dplyr::group_by(chromosome) %>%
    dplyr::summarize(
      begin = min(get(unit)),
      end = max(get(unit))
    ) %>%
    dplyr::ungroup()

  p +
    ggplot2::geom_segment(
      data = chromosomes,
      ggplot2::aes(x = begin + offset, xend = end - offset, y = y, yend = y),
      size = size, color = color, alpha = 1
    )

}
