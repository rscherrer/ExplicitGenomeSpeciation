#' Function to rename the columns of a data frame given a character vector
#'
#' @param df A data frame
#' @param newcols A character string with the new names of the columns
#'
#' @return A renamed data frame
#'
#' @examples
#'
#' df <- data.frame(a = 1, b = 1)
#' rename_str(df, c("aa", "bb"))
#'
#' @export

rename_str <- function(df, newcols) {

  if (length(newcols) != ncol(df)) {
    stop("Please provide the right number of column names")
  }
  colnames(df) <- newcols
  return (df)

}
