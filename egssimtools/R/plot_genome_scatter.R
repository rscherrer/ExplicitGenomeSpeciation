#' Genome scatter plot
#'
#' Flexible function to plot multiple aspects of the various loci against
#' each other. The data can be plotted for a given time point in a simulation,
#' or for no specific time point if the data are not timed (e.g. summary data).
#' The plot can be animated through time points. Edges of the gene network
#' can be displayed onto the layer of points representing the loci.
#' Supports additional aesthetics mapping and customization.
#'
#' @param data A data frame containing information on a per locus (and optionally per time point) basis
#' @param x,y Strings. Names of the variable on the horizontal and vertical axes, respectively
#' @param mode Either of "timepoint" to filter a specific time point, "summary" if there are no multiple time points, and "animation" if an animation through time points should be made
#' @param t Integer. The time point to filter
#' @param show_generation Logical. Whether to show the time point as a title
#' @param color,size Names of variables to map to those aesthetics (points).
#' @param by_trait Whether to facet by trait
#' @param edges A data frame containing information for each edge in the gene network, where connected loci are in columns names "from" and "to". Do not supply if the edges should not be shown.
#' @param edge_color,edge_size Names of variables to map to those aesthetics (edges).
#'
#' @return A ggplot, or an animation from the `gganimate` package
#'
#' @export

plot_genome_scatter <- function(

  data,
  x,
  y,
  mode = "timepoint",
  t = NULL,
  show_generation = FALSE,
  color = NULL,
  size = NULL,
  by_trait = FALSE,
  edges = NULL,
  edge_color = NULL,
  edge_size = NULL

) {

  if (mode == "summary") show_generation <- FALSE

  # Filter a single time point if needed
  if (mode == "timepoint") {

    if (is.null(t)) t <- dplyr::last(unique(data$time))
    data <- dplyr::filter(data, time == t)

  }

  # Set up
  p <- ggplot2::ggplot(data)

  # Plot network edges if needed
  if (!is.null(edges)) {

    if (by_trait) edges <- edges %>% dplyr::mutate(trait = paste("Trait", trait))

    if (mode == "animation") {

      # Edge coordinates for each time point in an animation
      edges <- data %>%
        dplyr::group_by(time) %>%
        tidyr::nest() %>%
        dplyr::mutate(edges = purrr::map(data, give_edges_coordinates, edges, x, y)) %>%
        dplyr::select(-data) %>%
        tidyr::unnest(cols = c(edges)) %>%
        dplyr::ungroup()

    } else {

      # Or a single set of edge coordinates
      edges <- give_edges_coordinates(data, edges, x, y)

    }

    aes_string <- "x = x0, xend = x1, y = y0, yend = y1"
    if (!is.null(edge_color)) aes_string <- paste0(aes_string, ", color = get(edge_color)")
    if (!is.null(edge_size)) aes_string <- paste0(aes_string, ", size = get(edge_size)")
    aes_string <- paste0("ggplot2::aes(", aes_string, ")")

    p <- p +
      ggplot2::geom_segment(
        data = edges, mapping = eval(rlang::parse_expr(aes_string)),
        size = 0.5, alpha = 0.5
      )

  }

  # Add the points after the edges
  aes_string <- "x = get(x), y = get(y)"
  if (!is.null(color)) aes_string <- paste0(aes_string, ", color = get(color)")
  if (!is.null(size)) aes_string <- paste0(aes_string, ", size = get(size)")
  aes_string <- paste0("ggplot2::aes(", aes_string, ")")
  p <- p + ggplot2::geom_point(mapping = eval(rlang::parse_expr(aes_string)))

  # Custom axis labels
  p <- p + xlab(x) + ylab(y)

  # Show the generation as title
  if (show_generation) {

    title <- paste(
      "Generation",
      ifelse(mode == "timepoint", t, "{round(frame_time, -2)}")
    )

    p <- p + ggplot2::labs(title = title)

  }

  # Facettize if needed
  if (by_trait) p <- ggsim::facettize(p, cols = "trait", prepend = "Trait ")

  # Animate if needed
  if (mode == "animation") p <- p + gganimate::transition_time(time)

  return(p)

}
