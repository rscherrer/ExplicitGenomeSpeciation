#' Set up parameter file
#'
#' Function to create a parameter file with the requested parameter values
#'
#' @param pars A named list or vector of parameter values. Numeric values will be coerced to character. Please provide compound parameters as strings with space-separated values. Names should be the names of the parameters.
#' @param template Optional file to add or change parameters into (the file is not overwritten).
#' @param saveto Optional file where to save the output.
#'
#' @return A vector of parameter entries, only if `saveto` is NULL.
#'
#' @export

set_param_file <- function(pars, template = NULL, saveto = NULL) {

  library(stringr)

  tmp <- character()
  if (!is.null(template)) tmp <- as.character(read.delim(template)[, 1])

  for (i in seq_along(pars)) {

    parname <- names(pars)[i]
    par <- pars[i]
    replacement <- paste(parname, par)

    # Search for the pattern in the file
    # If the pattern is found, replace all occurrences
    # Otherwise, add the new line to the end

    where <- grep(parname, tmp)
    if (length(where) > 0) {
      tmp[where] <- replacement
    } else{
      tmp <- c(tmp, replacement)
    }
  }

  if (is.null(saveto)) return (tmp)

  out <- file(saveto)
  writeLines(tmp, out)
  close(out)

}