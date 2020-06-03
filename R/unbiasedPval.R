#' Calculating Asymptotically Unbiased DP Two-Sided P-Value
#' @name pvalueTwoSide
#' @aliases pvaltwoside
#' @aliases pvalTwoSide
#'
#' @param Z A Binomial sample with Tulap noise
#' @param size The number of trials in Binomial distribution (parameter n in
#'   Binomial(n, \eqn{\theta}))
#' @param theta The success probability for each trial (parameter \eqn{\theta}
#'   in Binomial(n, \eqn{\theta}))
#' @param b Discrete Laplace noise parameters, obtained by \eqn{exp(-\epsilon)}
#' @param q The truncated quantiles
#'
#' @return A vector of asymptotically unbiased two-sided p-values
#' @seealso UMP one-sided p-values (\code{\link{pvalLeft}} and
#'   \code{\link{pvalRight}})
#' @references Awan, Jordan Alexander, and Aleksandra Slavkovic. 2020.
#'   "Differentially Private Inference for Binomial Data". Journal of Privacy
#'   and Confidentiality 10 (1). \url{https://doi.org/10.29012/jpc.725}.
#' @examples set.seed(2020)
#' sample <- rbinom(1, 10, 0.2) + rtulap(1, 0, 0.3, 0.05)
#' pvalTwoSide(sample, size = 10, theta = 0.5,b = 0.3, q = 0.05)
#'

#' @rdname pvalueTwoSide
#' @export
pvalTwoSide <- function(Z, size, theta, b, q){
  reps = base::length(Z)
  pval = base::rep(0,reps)
  values = base::seq(0,size)
  T = base::abs(Z-size*theta)
  return(pvalRight(Z = T+size*theta, size, theta, b, q) + 1-pvalRight(Z = size*theta-T, size, theta, b, q))
}
