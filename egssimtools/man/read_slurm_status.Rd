% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_slurm_status.R
\name{read_slurm_status}
\alias{read_slurm_status}
\title{Read the SLURM status of a simulation}
\usage{
read_slurm_status(root)
}
\arguments{
\item{root}{Path to the simulation folder}
}
\value{
A character
}
\description{
Checks the slurm output file
}
\examples{

\dontrun{

read_slurm_status("data/example_1") # should give an error because examples
# were not run using SLURM

}

}
