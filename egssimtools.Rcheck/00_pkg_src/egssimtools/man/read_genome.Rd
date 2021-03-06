% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_genome.R
\name{read_genome}
\alias{read_genome}
\title{Read locus-specific data through time}
\usage{
read_genome(
  folder,
  variables,
  parnames = NULL,
  combine = FALSE,
  as_numeric = NULL,
  architecture = FALSE,
  archfile = "architecture.txt",
  parfile = "paramlog.txt",
  nloci = NULL
)
}
\arguments{
\item{folder}{Path to the folder}

\item{variables}{What variables to read (\code{time} is included by default)}

\item{parnames, combine, as_numeric, architecture, archfile, parfile}{Parameters for \code{read_data}}

\item{nloci}{Number of loci (automatically guessed if unspecified)}
}
\value{
A data frame
}
\description{
Wrapper around \code{read_data} to read locus-wise data
}
\details{
The file \code{time.dat} must be present
}
\examples{

root <- system.file("extdata", "example_1", package = "egssimtools")

# Read Fst throughout the genome
read_genome(root, "genome_Fst")

# Read multiple metrics and attach architecture
read_genome(root, c("genome_Fst", "genome_Cst"), architecture = TRUE)

}
