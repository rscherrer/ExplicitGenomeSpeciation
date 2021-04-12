#' Plot gene network
#'
#' Show the gene interaction network as a graph.
#'
#' @param arch Genetic architecture, in `tbl_graph` format
#' @param data A data frame containing information on a per locus (and optionally per time point) basis
#' @param mode Either of "timepoint" to filter a specific time point, "summary" if there are no multiple time points, and "animation" if an animation through time points should be made
#' @param t Integer. The time point to filter
#' @param fill,size Names of the variables to map to these aesthetics
#' @param show_generation Logical. Whether to show the time point as a title.
#'
#' @return A ggplot
#'
#' @export

plot_genome_network <- function(

  arch,
  data,
  mode = "timepoint",
  t = NULL,
  fill = NULL,
  size = NULL,
  show_generation = FALSE

) {

  if (mode == "summary") show_generation <- FALSE

  # Filter a single time point if needed
  if (mode == "timepoint") {

    if (is.null(t)) t <- dplyr::last(unique(data$time))
    data <- dplyr::filter(data, time == t)

  }

  # Generate a layout for the graph
  layout <- igraph::layout_with_graphopt(
    arch, charge = 0.1, mass = 30, niter = 10000
  ) %>%
    as.data.frame() %>%
    dplyr::rename(x = V1, y = V2)

  # Add coordinates to the data
  data <- data %>%
    dplyr::group_by(time) %>%
    dplyr::mutate(x = layout$x, y = layout$y) %>%
    dplyr::ungroup()

  p <- ggraph::ggraph(arch, layout = layout) +
    ggraph::geom_edge_link(aes(color = trait), width = 0.2, alpha = 0.6)

  # Prepare aesthetic mapping
  aes_string <- "x = x, y = y"
  if (!is.null(fill)) aes_string <- paste0(aes_string, ", fill = get(fill)")
  if (!is.null(size)) aes_string <- paste0(aes_string, ", size = get(size)")
  aes_string <- paste0("ggplot2::aes(", aes_string, ")")

  # Add points
  p <- p + ggraph::geom_node_point(
      data = data,
      mapping = eval(rlang::parse_expr(aes_string)),
      shape = 21
    )

  # Show the generation as title
  if (show_generation) {

    title <- paste(
      "Generation",
      ifelse(mode == "timepoint", t, "{round(frame_time, -2)}")
    )

    p <- p + ggplot2::labs(title = title)

  }

  # Animate if needed
  if (mode == "animation") {

    p <- p + gganimate::transition_time(time)

  }

  return(p)

}
