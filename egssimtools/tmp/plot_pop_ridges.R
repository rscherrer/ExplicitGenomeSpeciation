#' Plot population distribution through time
#'
#' @param root Path to the simulation
#' @param y Variable to show
#' @param by Numbers of columns in which to split variable `y`. See `?read_data`.
#' @param j If variable `y` is splitted, which column to show?
#' @param times Optional vector of time points to show
#' @param facet_by Variable to facet by
#' @param cpp_numbering Parameter for `read_indiv_long`. Relevant only if
#' `facet_by` is specified.
#' @param ... Parameters to be passed to `ggridges::geom_density_ridges`
#'
#' @return A ggplot
#'
#' @examples
#'
#' \dontrun{
#'
#' root <- "data/example_1"
#'
#' # Plot the distribution of the ecological trait at three time points
#' plot_pop_ridges(
#'   root, "individual_trait", by = 3, j = 1, times = c(100, 200, 300)
#' )
#'
#' # Plot the distribution for the three traits
#' plot_pop_ridges(
#'   root, "individual_trait", by = 3, j = 1, times = c(100, 200, 300),
#'   facet_by = "trait"
#' )
#'
#' }
#'
#' @export

plot_pop_ridges <- function(
  root,
  y,
  by = 1,
  j = 1,
  times = NULL,
  facet_by = NULL,
  cpp_numbering = TRUE,
  ...
) {

  data <- read_pop(root, y, by = by)


  if (is.null(facet_by)) {

    ycol <- colnames(data)[grep(y, colnames(data))][j]
  } else {
    data <- read_indiv_long(
      root, y, by = by, names_to = facet_by, cpp_numbering = cpp_numbering
    )
    ycol <- "value"
  }

  if (is.null(times)) times <- unique(data$time)
  data <- data %>%
    dplyr::filter(time %in% times) %>%
    dplyr::mutate(time = factor(time))

  p <- ggplot2::ggplot(
    data,
    ggplot2::aes_string(x = ycol, y = "time")
  )

  if (!is.null(facet_by)) p <- p + ggplot2::aes(fill = factor(get(facet_by)))

  p <- p + ggridges:::geom_density_ridges(...) +
    ggplot2::theme_bw() +
    ggplot2::xlim(c(0, max(data[[ycol]]))) +
    ggplot2::theme(legend.position = "none")

  p %>% ggsim::facettize(cols = facet_by, prepend = paste0(facet_by, " "))

}
