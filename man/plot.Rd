\name{plot}
\alias{plot.tscausal}
\title{
Casual Impact Plots
}
\description{
Provides visual diagnostis for causal effects.
}
\usage{
\method{plot}{tscausal}(x, y = NULL, median_color = "black", 
median_type = 1, median_width = 2, gradient_color = "orange", 
interval_color = "red", interval_type = 2, interval_width = 2, 
ylim = NULL, ylab = "", n_original = NULL, ...)
}
\arguments{
\item{x}{
An object of class \dQuote{tscausal} generated from
calling \code{\link{tscausal}}.
}
\item{y}{
Not used.
}
\item{median_color}{
The color used for plotting the median value.
}
\item{median_type}{
The line type for the median.
}
\item{median_width}{
The width of the median line.
}
\item{gradient_color}{
The gradient color to use for the distribution.
}
\item{interval_color}{
The color of the quantile lines
}
\item{interval_type}{
The line type for the quantiles.
}
\item{interval_width}{
The width of the quantile lines.
}
\item{ylim}{
User specified limits for y.
}
\item{ylab}{
User specified label for y.
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
