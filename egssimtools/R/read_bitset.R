#' Read bitset
#'
#' Read a binary file bit-wise and return a vector of integers (0/1)
#' (e.g. whole genomes of the whole population)
#'
#' @param filename Path to the file to read
#'
#' @return A vector of integers
#'
#' @examples
#'
#' \dontrun{
#'
#' # Location of the simulation folder
#' root <- "data/example_1"
#'
#' f <- file.path(root, "time.dat")
#' read_bitset(f)
#'
#' }
#'
#' @export

read_bitset <- function(filename) {

  f <- file(filename, "rb")
  n <- file.size(filename)
  bytes <- readBin(f, raw(), n)
  close(f)
  return(as.integer(rawToBits(bytes)))

}
