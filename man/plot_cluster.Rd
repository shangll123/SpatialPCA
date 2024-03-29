% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpatialPCA_utilties.R
\name{plot_cluster}
\alias{plot_cluster}
\title{Visualize cluster labels on locations.}
\usage{
plot_cluster(
  location,
  clusterlabel,
  pointsize = 3,
  text_size = 15,
  title_in,
  color_in,
  legend = "none"
)
}
\arguments{
\item{location}{A n by k matrix of spot locations.}

\item{clusterlabel}{A vector of cluster labels for spots.}

\item{pointsize}{An integer, the point size of each spot.}

\item{title_in}{A character string, the title you want to display at the top of the figure.}

\item{color_in}{A vector of colors for each cluster.}

\item{legend}{A character string, the position of the figure legend. Select from "top", "bottom","left" or "right".}

\item{textsize}{An integer, the text size in the legend.}
}
\value{
A ggplot object.
}
\description{
This function visualizes cluster labels on locations.
}
