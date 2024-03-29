% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BINCO.r
\name{valley.score}
\alias{valley.score}
\title{6. valley point calculation #}
\usage{
valley.score(count.mix)
}
\description{
To calculate the valley score position, we first fit a smooth curve for the 
empirical density of selection frequencies using the R function
smooth.spline(df), where the df is determined such that the number of sign
changes of the smooth curve is 1. The valley point is at the position of 
the lowest of the curve.
}
