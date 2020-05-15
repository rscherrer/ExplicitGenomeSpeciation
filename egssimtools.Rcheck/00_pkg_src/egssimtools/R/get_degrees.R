#' Get locus degrees
#'
#' @param arch A genetic architecture
#'
#' @return A vector of degrees across loci
#'
#' @export

get_degrees <- function(arch) {

  library(tidyverse)

  degrees <- rep(0, length(arch$locations))
  connected <- table(do.call("c", arch$networks %>% map(~ do.call("c", .x[1:2]))))
  degrees[as.numeric(names(connected)) + 1] <- connected
  return(degrees)

}
