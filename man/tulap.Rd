% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tulap.R
\name{tulap}
\alias{tulap}
\alias{ptulap}
\alias{rtulap}
\title{The Truncated-Uniform-Laplace (Tulap) Distribution}
\usage{
rtulap(n, m = 0, b = 0, q = 0)

ptulap(t, m = 0, b = 0, q = 0)
}
\arguments{
\item{n}{number of observations}

\item{m}{vector of medians}

\item{b}{vector of Discrete Laplace noise parameters, obtained by
\eqn{exp(-\epsilon)}}

\item{q}{vector of truncated quantiles}

\item{t}{vector of quantiles}
}
\value{
\code{ptulap} gives the distribution function and \code{rtulap}
  generates random derviates.
}
\description{
Distribution function and random generation for the Tulap
  distribution.
}
\examples{
sample <- rtulap(1000, 1, 0.3, 0.05)
sample[1:10]
plot(density(sample))
ptulap(0.5, 0, 0.3, 0.05)

}
\references{
\code{ptulap} and \code{rtulap} are based on Awan, Jordan
  Alexander, and Aleksandra Slavkovic. 2020. "Differentially Private
  Inference for Binomial Data". Journal of Privacy and Confidentiality 10
  (1). \url{https://doi.org/10.29012/jpc.725}.
}
