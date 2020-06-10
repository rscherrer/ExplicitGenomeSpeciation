#' String to numeric vector
#'
#' @param string A string containing numbers
#' @param sep The separator between numbers
#'
#' @return A numeric vector

str2vec <- function(string, sep = ' ') {

  as.numeric(unlist(strsplit(as.character(string), sep)))

}
