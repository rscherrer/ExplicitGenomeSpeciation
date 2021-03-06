% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_genome_scatter.R
\name{plot_genome_scatter}
\alias{plot_genome_scatter}
\title{Genome scatter plot}
\usage{
plot_genome_scatter(
  data,
  x,
  y,
  mode = "timepoint",
  t = NULL,
  show_generation = FALSE,
  color = NULL,
  size = NULL,
  by_trait = FALSE,
  edges = NULL,
  edge_color = NULL,
  edge_size = NULL
)
}
\arguments{
\item{data}{A data frame containing information on a per locus (and optionally per time point) basis}

\item{x, y}{Strings. Names of the variable on the horizontal and vertical axes, respectively}

\item{mode}{Either of "timepoint" to filter a specific time point, "summary" if there are no multiple time points, and "animation" if an animation through time points should be made}

\item{t}{Integer. The time point to filter}

\item{show_generation}{Logical. Whether to show the time point as a title}

\item{color, size}{Names of variables to map to those aesthetics (points).}

\item{by_trait}{Whether to facet by trait}

\item{edges}{A data frame containing information for each edge in the gene network, where connected loci are in columns names "from" and "to". Do not supply if the edges should not be shown.}

\item{edge_color, edge_size}{Names of variables to map to those aesthetics (edges).}
}
\value{
A ggplot, or an animation from the \code{gganimate} package
}
\description{
Flexible function to plot multiple aspects of the various loci against
each other. The data can be plotted for a given time point in a simulation,
or for no specific time point if the data are not timed (e.g. summary data).
The plot can be animated through time points. Edges of the gene network
can be displayed onto the layer of points representing the loci.
Supports additional aesthetics mapping and customization.
}
