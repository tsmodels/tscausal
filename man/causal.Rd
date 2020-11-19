\name{tscausal}
\alias{tscausal}
\alias{tscausal.tsmodel.distribution}
\title{
Time Series Based Causal Analysis
}
\description{
Performs a time series causal analysis using a simulated predictive distribution
for the counterfactual.
}
\usage{
tscausal(object, \dots)
\method{tscausal}{tsmodel.distribution}(object, actual, fitted = NULL, 
alpha = 0.05, \dots)
}
\arguments{
\item{object}{
An object which inherits class \dQuote{tsmodel.distribution} representing
the forecast distribution in the post intervention period.
}
\item{actual}{
An xts vector of the actual data series which includes both the pre
and post intervention data.
}
\item{fitted}{
An optional object for the in-sample fitted values which can be either
an xts vector or an object of class \dQuote{tsmodel.distribution}.
}
\item{alpha}{
The coverage representing the 1-alpha confidence level.
}
\item{\ldots}{
Any additional arguments passed to custom classes.
}
}
\details{
The routine calculates the point wise differences between the actual and
counterfactual (distribution) to determine the distribution of the lift.
}
\value{
An object of class \dQuote{tscausal} which can be passed to the \code{print},
\code{report} or \code{plot} functions.
}
\references{
Brodersen, Kay H and Gallusser, Fabian and Koehler, Jim and Remy, Nicolas and
Scott, Steven L and others (2016). Inferring causal impact using Bayesian structural
time-series models. The Annals of Applied Statistics, 9 (1), 247-274.
}
\author{
Alexios Galanos with some supporting code borrowed from the CausalImpact
package of Scott.
}
