#' Make facet labels
#'
#' Make facet labels of the form "variable = value".
#'
#' @param data The dataset
#' @param varname The variable used for facetting
#'
#' @export

make_facet_labels <- function(data, varname) {

  labs <- unique(unlist(data[, varname]))
  labnames <- as.character(labs)
  labs <- paste(varname, labs, sep = " = ")
  names(labs) <- labnames
  return (labs)

}
