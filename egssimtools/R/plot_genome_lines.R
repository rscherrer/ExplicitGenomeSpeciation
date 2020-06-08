#' Plot locus variables through time
#'
#' @param root Path to the simulation
#' @param y Variable to plot
#' @param span Optional span parameter for LOESS smoothing (no smoothing if 0)
#' @param facet_by Variable to facet by
#' @param archfile Optional architecture file name
#'
#' @return A ggplot
#'
#' @examples
#'
#' root <- system.file("extdata", "example_1", package = "egssimtools")
#' plot_genome_lines(root, "genome_Fst")
#' plot_genome_lines(root, "genome_Fst", facet_by = "trait")
#'
#' @export

plot_genome_lines <- function(
  root,
  y,
  span = 0.2,
  facet_by = NULL,
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
  )
  if (!is.null(facet_by)) {
    p <- p +
      ggplot2::aes(color = factor(trait)) +
      ggplot2::scale_color_manual(values = colors)
  }
  p <- p +
    ggplot2::labs(x = "time", y = y) +
    ggplot2::theme(legend.position = "none")
  if (!is.null(facet_by)) p <- p + ggplot2::labs(color = facet_by)

  p <- p + ggplot2::ylim(c(0, max(data[[y]])))

  p %>% ggsim::facettize(rows = facet_by, prepend = paste0(facet_by, " "))

}
