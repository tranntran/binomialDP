#' @title Calculating UMP One-Sided p-values
#' @name pvalueOneSide
#' @aliases pval
#' @aliases pvalleft
#' @aliases pvalright
#'
#' @description Calculating UMP one-sided p-values for binomial data under
#'   \eqn{(\epsilon, \delta)}-DP
#' @param Z A Binomial sample with Tulap noise
#' @param size The number of trials in Binomial distribution (parameter n in
#'   Binomial(n, \eqn{\theta}))
#' @param theta The success probability for each trial (parameter \eqn{\theta}
#'   in Binomial(n, \eqn{\theta}))
#' @param b Discrete Laplace noise parameters, obtained by \eqn{exp(-\epsilon)}
#' @param q The truncated quantiles
#'
#' @return A vector of one-sided p-values
#' @seealso Asymptotically unbiased DP two-sided p-values (\code{\link{pvalTwoSide}})
#' @references Awan, Jordan Alexander, and Aleksandra Slavkovic. 2020.
#'   "Differentially Private Inference for Binomial Data". Journal of Privacy
#'   and Confidentiality 10 (1). \url{https://doi.org/10.29012/jpc.725}.
#' @examples
#' set.seed(2020)
#' sample <- rbinom(1, 10, 0.2) + rtulap(1, 0, 0.3, 0.05)
#' pvalLeft(sample, size = 10, theta = 0.5,b = 0.3, q = 0.05)
#' pvalRight(sample, size = 10, theta = 0.5,b = 0.3, q = 0.05)
NULL


#' @rdname pvalueOneSide
#' @export
pvalRight <- function(Z, size, theta, b, q){
  reps = base::length(Z)
  pval = base::rep(0,reps)
  values = base::seq(0,size)

  B = stats::dbinom(values, size = size, prob = theta)

  for(r in 1:reps){
    #F = ptulap(t = Z[r]-values, m = 0, b, q)
    F = ptulap(t = values-Z[r], m = 0, b, q)
    pval[r]=t(F)%*%B
  }
  return(pval)
}

#' @rdname pvalueOneSide
#' @export
pvalLeft <- function(Z, size, theta, b, q){
  reps = base::length(Z)
  pval = base::rep(0,reps)
  values = base::seq(0,size)

  B = stats::dbinom(values, size = size, prob = theta)

  for(r in 1:reps){
    F = 1 - ptulap(t = values-Z[r], m = 0, b, q)
    #F = ptulap(t = values-Z[r], m = 0, b, q)
    pval[r]=t(F)%*%B
  }
  return(pval)
}
