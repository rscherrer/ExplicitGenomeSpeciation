#' Read an architecture file
#'
#' @param filename Path to the architecture file
#'
#' @return A list with all the fields of a genetic architecture
#'
#' @examples
#'
#' \dontrun{
#'
#' # Location of the simulation folder
#' root <- "data/example_1"
#'
#' f <- file.path(root, "architecture.txt")
#' read_arch_file(f)
#'
#' }
#'
#' @export

read_arch_file <- function(filename) {

  archfile <- readLines(filename)

  begin <- grep("^Architecture:$", archfile)
  if (length(begin) == 0) stop(
    paste0(
      "Cannot find the beginning of the genetic architecture in file ",
      filename
    )
  )
  begin <- begin + 1
  archfile <- archfile[begin:length(archfile)]

  fields <- c("chromosomes", "traits", "locations", "effects", "dominances")
  arch <- lapply(fields, function(x) {
    str2vec(archfile[which(archfile == x)[1] + 1])
  })
  names(arch) <- fields

  ii <- grep("network", archfile)
  networks <- lapply(ii, function(i) list (
    edges0 = egssimtools:::str2vec(archfile[i + 2]),
    edges1 = egssimtools:::str2vec(archfile[i + 4]),
    weights = egssimtools:::str2vec(archfile[i + 6])
  ))
  names(networks) <- c("ecotrait", "matepref", "neutral")
  arch$networks <- networks

  return(arch)

}
