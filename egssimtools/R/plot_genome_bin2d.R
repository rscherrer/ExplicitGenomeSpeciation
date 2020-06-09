#' Plot genome distribution through time
#'
#' @param root Path to the simulation folder
#' @param y Variable to plot the density for
#' @param x Optional variable for the horizontal axis. Defaults to "time".
#' @param t What time point? Defaults to last time point, only if `x` is not "time"
#' @param facet_by Optional variable to facet by
#' @param archfile Optional architecture file name
#' @param ... Parameters to be passed to `geom_path`
#'
#' @return A ggplot
#'
#' @examples
#'
#' root <- system.file("extdata", "example_1", package = "egssimtools")
#' plot_genome_bin2d(root, "genome_Fst")
#'
#' # Plot Fst versus Cst
#' plot_genome_bin2d(root, "genome_Fst", x = "genome_Cst")
#'
#' # Can facet by some variable
#' plot_genome_bin2d(root, "genome_Fst", facet_by = "trait")
#'
#' @export

plot_genome_bin2d <- function(
  root,
  y,
  x = "time",
  t = NULL,
  facet_by = NULL,
  archfile = "architecture.txt",
  ...
) {

  to_read <- y
  if (x != "time" & x != y) to_read <- c(x, y)
  data <- read_loci(root, to_read, architecture = TRUE, archfile = archfile)

  if (x != "time") {
    if (is.null(t)) t <- dplyr::last(data$time)
    data <- data %>% dplyr::filter(time == t)
  }

  p <- ggplot2::ggplot(data, ggplot2::aes_string(x = x, y = y)) +
    ggplot2::geom_bin2d(...) +
    ggplot2::theme_bw()

  p %>% ggsim::facettize(cols = facet_by, prepend = paste0(facet_by, " "))

}
