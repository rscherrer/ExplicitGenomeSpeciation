#' Plot genome distribution through time
#'
#' @param root Path to the simulation
#' @param y Variable to show
#' @param times Optional vector of time points to show
#' @param facet_by Variable to facet by
#' @param archfile Optional architecture file name
#'
#' @return A ggplot
#'
#' @examples
#'
#' root <- system.file("extdata", "example_1", package = "egssimtools")
#' plot_genome_ridges(root, "genome_Fst", times = c(100, 200, 300))
#' plot_genome_ridges(
#'  root, "genome_Fst", times = c(100, 200, 300), facet_by = "trait"
#' )
#'
#' @export

plot_genome_ridges <- function(
  root,
  y,
  times = NULL,
  facet_by = NULL,
  archfile = "architecture.txt"
) {

  data <- read_loci(root, y, architecture = TRUE, archfile = archfile)

  if (is.null(times)) times <- unique(data$time)
  data <- data %>% dplyr::filter(time %in% times)

  colors <- c("forestgreen", "goldenrod", "lightgrey")

  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = get(y), y = factor(time))
  )
  if (!is.null(facet_by)) p <- p + ggplot2::aes(fill = factor(get(facet_by)))
  p <- p + ggridges:::geom_density_ridges() +
    ggplot2::theme_bw() +
    ggplot2::labs(x = y, y = "time") +
    ggplot2::xlim(c(0, max(data[[y]]))) +
    ggplot2::theme(legend.position = "none")
  if (!is.null(facet_by)) {
    p <- p +
      ggplot2::labs(fill = facet_by) +
      ggplot2::scale_fill_manual(values = colors)
  }

  p %>% ggsim::facettize(cols = facet_by, prepend = paste0(facet_by, " "))

}
