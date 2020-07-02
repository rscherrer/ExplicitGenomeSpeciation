#' Read bitset
#'
#' Read a binary file bit-wise and return a vector of integers (0/1)
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

  # Read raw bytes from file, convert them into bits then integers (0/1)
  x <- readBin(file(filename, "rb"), raw(), n = file.size(filename))
  return(as.integer(rawToBits(x)))

}
