#' Genome-wide distributions through time
#'
#' @param root Path to the simulation
#' @param y Variable to show
#' @param times Optional vector of time points to show
#'
#' @return A ggplot
#'
#' @export

plot_genome_ridges <- function(root, y, times = NULL) {

  library(ggplot2)
  library(ggridges)
  library(ggsim)

  time <- read_data(root, "time")
  X <- read_data(root, y)
  nloci <- nrow(X) / nrow(time)

  data <- read_data(root, c("time", y), dupl = c(nloci, 1), architecture = TRUE)

  if (is.null(times)) times <- unique(data$time)
  data <- data %>% filter(time %in% times)

  p <- ggplot(data, aes(x = get(y), y = factor(time), fill = factor(trait))) +
    geom_density_ridges() +
    theme_bw() +
    labs(x = y, y = "Time", fill = "Trait") +
    scale_fill_manual(values = c("forestgreen", "goldenrod", "lightgrey")) +
    xlim(c(0, max(data[[y]]))) +
    theme(legend.position = "none")

  p %>% facettize(cols = "trait", prepend = "Trait ")

}
