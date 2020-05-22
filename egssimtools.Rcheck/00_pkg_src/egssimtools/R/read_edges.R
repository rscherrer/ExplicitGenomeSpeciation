#' Read edge-specific data
#'
#' @param filename Path to the folder
#' @param variable What to read
#'
#' @export

read_edges <- function(folder, variable) {

  t <- read_time(folder)
  data <- read_data(folder, variable)
  nedges <- length(data) / length(t)
  timepoints <- rep(t, each = nedges)
  data <- split(data, timepoints)

  # But how many edges per trait?
  nedges_per_trait <- as.numeric(read_parameters(folder, "nedges")$nedges)
  nedges_per_trait
  traits <- mrep(0:2, nedges_per_trait)

  lapply(data, split, traits) # returns a list of generation-lists of trait-vectors of edge properties

}
