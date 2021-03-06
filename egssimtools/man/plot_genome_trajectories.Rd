% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_genome_trajectories.R
\name{plot_genome_trajectories}
\alias{plot_genome_trajectories}
\title{Plot genomic trajectories}
\usage{
plot_genome_trajectories(data, y, color = NULL, alpha = 0.5, by_trait = FALSE)
}
\arguments{
\item{data}{A locus-wise data frame with columns "time", "locus" and \code{y} at least}

\item{y}{Name of the variable to plot}

\item{color}{Name of the variable to map to that aesthetic.}

\item{alpha}{Transparency value}

\item{by_trait}{Whether to facet by trait}
}
\value{
A ggplot
}
\description{
Trajectories of loci through time along a certain variable.
}
