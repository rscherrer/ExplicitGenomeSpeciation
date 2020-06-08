#' Genome-wide distributions through time
#'
#' @param root Path to the simulation
#' @param y Variable to show
#' @param times Optional vector of time points to show
#' @param archfile Optional architecture file name
#'
#' @return A ggplot
#'
#' @examples
#'
#' root <- system.file("extdata", "example_1", package = "egssimtools")
#' plot_genome_ridges(root, "genome_Fst")
#'
#' @export

plot_genome_ridges <- function(
  root,
  y,
  times = NULL,
  archfile = "architecture.txt"
) {

  data <- read_loci(root, y, architecture = TRUE, archfile = archfile)

  if (is.null(times)) times <- unique(data$time)
  data <- data %>% dplyr::filter(time %in% times)

  colors <- c("forestgreen", "goldenrod", "lightgrey")

  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = get(y), y = factor(time), fill = factor(trait))
  ) +
    ggridges:::geom_density_ridges() +
    ggplot2::theme_bw() +
    ggplot2::labs(x = y, y = "Time", fill = "Trait") +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::xlim(c(0, max(data[[y]]))) +
    ggplot2::theme(legend.position = "none")

  p %>% ggsim::facettize(cols = "trait", prepend = "Trait ")

}
