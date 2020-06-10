#' Read locus-specific genome architecture
#'
#' Make a table from locus-specific genetic architecture details
#'
#' @param folder Path to the simulation
#' @param filename Name of the architecture file
#'
#' @return A data frame with, for each locus, its location, trait, effect, dominance, chromosome and degree
#'
#' @examples
#'
#' root <- system.file("extdata", "example_1", package = "egssimtools")
#' read_arch_genome(root)
#'
#' @export

read_arch_genome <- function(folder, filename = "architecture.txt") {

  arch <- read_arch(folder, filename)

  data.frame(
    locus = seq(arch$location),
    location = arch$locations,
    trait = factor(arch$traits),
    effect = arch$effects,
    dominance = arch$dominances,
    chromosome = get_chromosomes(arch),
    degree = get_degrees(arch)
  )

}
