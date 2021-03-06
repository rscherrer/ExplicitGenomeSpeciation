% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_chromosomes.R
\name{get_chromosomes}
\alias{get_chromosomes}
\title{Get locus chromosomal locations}
\usage{
get_chromosomes(arch)
}
\arguments{
\item{arch}{A genetic architecture}
}
\value{
A vector of chromosome indices for each locus
}
\description{
Get locus chromosomal locations
}
\examples{

\dontrun{

# Location of the simulation folder
root <- "data/example_1"

arch <- read_arch(root)
get_chromosomes(arch)

}

}
