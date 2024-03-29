% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BINCO.r
\name{BINCO.plot}
\alias{BINCO.plot}
\title{check the null estimate by plot}
\usage{
BINCO.plot(out.BINCO, show.range = NULL, show.cutoff = TRUE, ...)
}
\description{
this function provides an overlay plot of the empirical mixture distribution of selection frequencies and the null distribution estimated by BINCO
out.BINCO is the outcome retured from BINCO function.
show.range is a vector of length=2, which defines the display range of the overlay plot. We recommand that the lower bound of this range should be greater than the lower bound of the fitting range used by BINCO.
show.cutoff determines whether a vertical line will be drawn at the valley point value calculated by BINCO
}
