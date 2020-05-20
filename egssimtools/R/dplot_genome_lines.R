#' Plot locus-specific lines through time
#'
#' @param roooooooot Path to the simulation
#' @param y Variable to plot
#' @param span Optional span parameter for LOESS smoothing (no smoothing if 0)
#'
#' @return A ggplot
#'
#' @export

dplot_genome_lines <- function(root, y, span = 0.2) {

  library(ggplot2)
  library(ggridges)
  library(ggsim)

  time <- read_data(root, "time")
  X <- read_data(root, y)
  nloci <- nrow(X) / nrow(time)

  data <- read_data(root, c("time", y), dupl = c(nloci, 1), architecture = TRUE)

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
    group_by(locus) %>%
    nest() %>%
    mutate(smooth = map(data, smoothen)) %>%
    unnest(cols = c(data, smooth)) %>%
    ungroup()

  p <- gglineplot(data, x = "time", y = "smooth", line = "locus", alpha = 0.5) +
    aes(color = factor(trait)) +
    scale_color_manual(values = c("forestgreen", "goldenrod", "lightgrey")) +
    labs(x = "Time (generations)", y = y, color = "Trait") +
    theme(legend.position = "none")

  p <- p + ylim(c(0, max(data[[y]])))

  p %>% facettize(rows = "trait", prepend = "Trait ")

}
