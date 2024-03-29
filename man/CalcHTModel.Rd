% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbtm_functions.R
\name{calcHTModel}
\alias{calcHTModel}
\title{A Function to calculate the Hydrotime model parameters.}
\usage{
calcHTModel(
  data,
  germ.wp = "GermWP",
  cum.time = "CumTime",
  cum.frac = "CumFraction",
  max.cum.frac = 1,
  plot = TRUE
)
}
\arguments{
\item{data}{time course and cumulative dataset to be used in the hydrotime model. The original dataframe template should be used or column names should be modified similarly to the template. A column with time in hours CumTime) + a column with cumulative fractions (CumFraction) and the experiment water potential (Germ.wp) are required. Filter the dataframe to only have treatments with water potential at the same temperature level.}

\item{germ.wp}{Column containing germination water potential.}

\item{cum.time}{Column containing cumulative elapsed time.}

\item{cum.frac}{Column containing cumulative fraction germinated.}

\item{max.cum.frac}{Sets the ceiling cumulative fraction for the model when treatment at optimal condition displays a lower maximum cumulative fraction. Use it on your own discretion.}

\item{plot}{Should the model results be plotted?}
}
\description{
This function calculates the hydrotime constant (HT), the median water potential base (Psib50) and the standard deviation (sigma).
}
\examples{
"foo"
}
\keyword{hydrotime}
\keyword{model}
\keyword{parameters}
