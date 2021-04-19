#' Get locus degrees
#'
#' @param arch A genetic architecture
#'
#' @return A vector of degrees across loci
#'
#' @examples
#'
#' \dontrun{
#'
#' # Location of the simulation folder
#' root <- "data/example_1"
#'
#' arch <- read_arch(root)
#' get_degrees(arch)
#'
#' }
#'
#' @export

get_degrees <- function(arch) {

  degrees <- rep(0, length(arch$locations))
  connected <- table(
    do.call("c", purrr::map(arch$networks, ~ do.call("c", .x[1:2])))
  )
  degrees[as.numeric(names(connected)) + 1] <- connected
  # numbering starts at zero in C++
  return(degrees)

}
