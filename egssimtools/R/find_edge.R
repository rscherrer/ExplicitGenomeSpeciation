#' Find an edge in a gene network
#'
#' Returns a logical vector indicating which is the target edge, from a series of edges in a genetic architecture
#'
#' @param i,j The indices of the loci connected by the edge
#' @param arch A genetic architecture dataset, on a per-edge basis (see e.g. `?read_arch_network`)
#'
#' @return A logical vector
#'
#' @export

find_edge <- function(i, j, arch) {

  with(arch, from == i & to == j | from == j & to == i)

}
