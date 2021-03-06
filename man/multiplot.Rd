% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reportfunctions.R
\name{multiplot}
\alias{multiplot}
\title{Multiple plot function}
\usage{
multiplot(..., plotlist = NULL, cols = 1, layout = NULL)
}
\arguments{
\item{cols:}{Number of columns in layout}

\item{layout:}{A matrix specifying the layout. If present, 'cols' is ignored.}

\item{plotlist:}{vector of plots}
}
\description{
ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
then plot 1 will go in the upper left, 2 will go in the upper right, and
3 will go all the way across the bottom.
}
\author{
Stefan Bartels, \email{email@biobits.eu}
}
