% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_arch.R
\name{read_arch}
\alias{read_arch}
\title{Read architecture}
\usage{
read_arch(folder, filename = "architecture.txt")
}
\arguments{
\item{folder}{Path to the folder}

\item{filename}{Optional file name}
}
\value{
A list with all the fields of a genetic architecture
}
\description{
Read architecture
}
\examples{

root <- system.file("extdata", "example_1", package = "egssimtools")
read_arch(root)

}
