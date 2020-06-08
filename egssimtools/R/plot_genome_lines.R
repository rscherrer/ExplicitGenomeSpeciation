#' Plot locus-specific lines through time
#'
#' @param root Path to the simulation
#' @param y Variable to plot
#' @param span Optional span parameter for LOESS smoothing (no smoothing if 0)
#' @param archfile Optional architecture file name
#'
#' @return A ggplot
#'
#' @export

plot_genome_lines <- function(
  root,
  y,
  span = 0.2,
  archfile = "architecture.txt"
) {

  data <- read_loci(root, y, architecture = TRUE, archfile = archfile)

  if (span > 0) {
    smoothen <- function(data) {
      loess(
        as.formula(paste(y, "~ time")), degree = 1, span = span, data = data
      )$fitted
    }
  } else {
    smoothen <- function(data) data[[y]]
  }

  data <- data %>%
    dplyr::group_by(locus) %>%
    tidyr::nest() %>%
    dplyr::mutate(smooth = purrr::map(data, smoothen)) %>%
    tidyr::unnest(cols = c(data, smooth)) %>%
    dplyr::ungroup()

  colors <- c("forestgreen", "goldenrod", "lightgrey")

  p <- ggsim::gglineplot(
    data,
    x = "time",
    y = "smooth",
    line = "locus",
    alpha = 0.5
  ) +
    ggplot2::aes(color = factor(trait)) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(x = "Time (generations)", y = y, color = "Trait") +
    ggplot2::theme(legend.position = "none")

  p <- p + ggplot2::ylim(c(0, max(data[[y]])))

  p %>% ggsim::facettize(rows = "trait", prepend = "Trait ")

}
