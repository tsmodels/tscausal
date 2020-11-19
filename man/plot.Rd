\name{plot}
\alias{plot.tscausal}
\title{
Casual Impact Plots
}
\description{
Provides visual diagnostis for causal effects.
}
\usage{
\method{plot}{tscausal}(x, y = NULL, n_original = 100, ...)
}
\arguments{
\item{x}{
An object of class \dQuote{tscausal} generated from
calling \code{\link{tscausal}}.
}
\item{y}{
Not used.
}
\item{n_original}{
Number of points from the end of the intervention period
to include in the plotting. For large time series it makes
sense to not plot the entire period pre-intervention as
this will visually obscure the post-intervention effects.
}
\item{\dots}{
Additional parameters passed to the \dQuote{plot.tsmodel.distribution}
functions in the \dQuote{tsmethods} package.
}
}
\author{
Alexios Galanos
}
