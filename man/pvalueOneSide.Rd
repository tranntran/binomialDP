% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pvalue.R
\name{pvalueOneSide}
\alias{pvalueOneSide}
\alias{pval}
\alias{pvalleft}
\alias{pvalright}
\alias{pvalRight}
\alias{pvalLeft}
\title{Calculating UMP One-Sided p-values}
\usage{
pvalRight(Z, size, theta, b, q)

pvalLeft(Z, size, theta, b, q)
}
\arguments{
\item{Z}{A Binomial sample with Tulap noise}

\item{size}{The number of trials in Binomial distribution (parameter n in
Binomial(n, \eqn{\theta}))}

\item{theta}{The success probability for each trial (parameter \eqn{\theta}
in Binomial(n, \eqn{\theta}))}

\item{b}{Discrete Laplace noise parameters, obtained by \eqn{exp(-\epsilon)}}

\item{q}{The truncated quantiles}
}
\value{
A vector of one-sided p-values
}
\description{
Calculating UMP one-sided p-values for binomial data under
  \eqn{(\epsilon, \delta)}-DP
}
\examples{
set.seed(2020)
sample <- rbinom(1, 10, 0.2) + rtulap(1, 0, 0.3, 0.05)
pvalLeft(sample, size = 10, theta = 0.5,b = 0.3, q = 0.05)
pvalRight(sample, size = 10, theta = 0.5,b = 0.3, q = 0.05)
}
\references{
Awan, Jordan Alexander, and Aleksandra Slavkovic. 2020.
  "Differentially Private Inference for Binomial Data". Journal of Privacy
  and Confidentiality 10 (1). \url{https://doi.org/10.29012/jpc.725}.
}
\seealso{
Asymptotically unbiased DP two-sided p-values (\code{\link{pvalTwoSide}})
}
