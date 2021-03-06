% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_genome_network.R
\name{plot_genome_network}
\alias{plot_genome_network}
\title{Plot gene network}
\usage{
plot_genome_network(
  arch,
  data,
  mode = "timepoint",
  t = NULL,
  fill = NULL,
  size = NULL,
  show_generation = FALSE
)
}
\arguments{
\item{arch}{Genetic architecture, in \code{tbl_graph} format}

\item{data}{A data frame containing information on a per locus (and optionally per time point) basis}

\item{mode}{Either of "timepoint" to filter a specific time point, "summary" if there are no multiple time points, and "animation" if an animation through time points should be made}

\item{t}{Integer. The time point to filter}

\item{fill, size}{Names of the variables to map to these aesthetics}

\item{show_generation}{Logical. Whether to show the time point as a title.}
}
\value{
A ggplot
}
\description{
Show the gene interaction network as a graph.
}
