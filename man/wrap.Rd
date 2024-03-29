% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/activity_code.r
\name{wrap}
\alias{wrap}
\title{Wraps data on a given range.}
\usage{
wrap(x, bounds = c(0, 2 * pi))
}
\arguments{
\item{x}{A vector of numeric data.}

\item{bounds}{The range within which to wrap \code{x} values}
}
\value{
A vector of numeric values within the limits defined by \code{bounds}
}
\description{
Input data outside the given bounds (default radian [0, 2*pi]) are wrapped to appear within the range.
}
\details{
As an example of wrapping, on bounds [0, 1], a value of 1.2 will be converted to 0.2, while a value of -0.2 will be converted to 0.8.
}
\examples{
data(BCItime)
adjtime <- BCItime$time + 1/24
summary(adjtime)
adjtime <- wrap(adjtime, c(0,1))
summary(adjtime)
}
