#' Read binary file
#'
#' @param filename Path to the file
#' @param digits Encoding
#'
#' @return A numeric vector
#'
#' @export

read_binary <- function(filename, digits = 8) {

  # Number of values to read
  n <- file.info(filename)$size / digits

  file <- file(filename, "rb")
  data <- readBin(file, numeric(), n)
  close(file)

  return (data)

}
