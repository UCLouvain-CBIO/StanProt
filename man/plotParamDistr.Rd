% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ProtModelPlots.R
\name{plotParamDistr}
\alias{plotParamDistr}
\title{Plots and compare the distributions for a given parameter}
\usage{
plotParamDistr(
  x1,
  x2 = NULL,
  x3 = NULL,
  paramName,
  xlim = NULL,
  ylim = NULL,
  keepBinWidth = TRUE,
  vline_xintercept = NULL
)
}
\arguments{
\item{x1}{A \code{ProtModel} object which should have run already}

\item{x2}{An optional \code{ProtModel} object which should have run already}

\item{paramName}{A \code{character} containing the name of the parameter to plot}
}
\value{
a ggplot graph showing both distributions for the given parameter
}
\description{
Plot either one or two distributions for a given parameter.
}
\author{
Philippe Hauchamps
}
