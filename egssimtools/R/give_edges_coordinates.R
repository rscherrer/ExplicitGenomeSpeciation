#' Attach locus coordinates to edges
#'
#' Add x and y corrdinates to each locus of each edge in an edge-wise
#' dataset.
#'
#' @param data A locus-wise data frame. It needs at least a column "locus", and columns named after arguments `x` and `y`.
#' @param edges An edge-wise data frame, with at least a column "from"and a column "to".
#' @param x,y Names of the columns of `data` that contain locus coordinates.
#'
#' @return An edge-wise data frame
#'
#' @export

give_edges_coordinates <- function(data, edges, x, y) {

  edges %>%
    dplyr::mutate(
      x0 = purrr::map_dbl(from, ~ with(data, get(x)[locus == .x])),
      y0 = purrr::map_dbl(from, ~ with(data, get(y)[locus == .x])),
      x1 = purrr::map_dbl(to, ~ with(data, get(x)[locus == .x])),
      y1 = purrr::map_dbl(to, ~ with(data, get(y)[locus == .x]))
    )

}
