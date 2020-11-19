\name{tsreport}
\alias{tsreport.tscausal}
\alias{tsreport}
\title{
Report Causal Impact Analysis
}
\description{
Provides a verbose analysis of the causal impact.
}
\usage{
\method{tsreport}{tscausal}(object, digits = 4, doc_template = NULL, 
type = c("screen","pdf","doc","html"), output_dir = "/", 
args = list(name = "Causal Analysis", frequency = NULL, model = NULL), ...)
}
\arguments{
\item{object}{
An object of class \dQuote{tscausal} generated from
calling \code{\link{tscausal}}.
}
\item{digits}{
Number of digits to print.
}
\item{doc_template}{
Filename (full path) of a docx which can be used as a template for the word document
option.
}
\item{type}{
The type of output.
}
\item{output_dir}{
A directory where the output is written and the causal object saved.
}
\item{args}{
An additional list of arguments used in the generation of the report.
}
\item{\dots}{
Not currently used.
}
}
\value{
Text output
}
\references{
Brodersen, Kay H and Gallusser, Fabian and Koehler, Jim and Remy, Nicolas and
Scott, Steven L and others (2016). Inferring causal impact using Bayesian structural
time-series models. The Annals of Applied Statistics, 9 (1), 247-274.
}
\author{
Alexios Galanos with supporting code borrowed from the CausalImpact
package of Scott.
}
