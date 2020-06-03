#' @title Finding One-Sided Confidence Intervals
#' @name CIOneSide
#' @aliases cioneside
#' @aliases cilower
#' @aliases ciupper
#'
#' @param alpha The confidence level of the returned confidence intervals
#' @param Z A Binomial sample with Tulap noise
#' @param size The number of trials in Binomial distribution (parameter n in
#'   Binomial(n, \eqn{\theta}))
#' @param b Discrete Laplace noise parameters, obtained by \eqn{exp(-\epsilon)}
#' @param q The truncated quantiles
#'
#' @return `CILower` returns a lower bound while `CIUpper` returns an upper bound.
#' @references Awan, Jordan Alexander, and Aleksandra Slavkovic. 2020.
#'   "Differentially Private Inference for Binomial Data". Journal of Privacy
#'   and Confidentiality 10 (1). \url{https://doi.org/10.29012/jpc.725}.
#' @seealso Finding asymptotically unbiased two-sided confidence intervals
#'   \code{\link{CITwoSide}}
#'
#' @examples
#' CILower(0.05, 7, 10, exp(-1), 0.05) #confidence interval (0.384, 1)
#' CIUpper(0.05, 2, 10, exp(-1), 0.05) #confidence interval (0, 0.517)
NULL

#' @rdname CIOneSide
#' @export
CILower <- function(alpha, Z, size, b, q){
  CIobj = function(theta, alpha, Z, size, b, q){
    return((pvalRight(Z = Z, size = size, theta = theta, b = b, q = q) - alpha)^2)
  }
  L = stats::optim(par = .5, fn = CIobj,
                   alpha = alpha, Z = Z, size = size, b = b, q = q,
                   method = "Brent", lower = 0, upper = 1)
  return(L$par)
}


#' @rdname CIOneSide
#' @export
CIUpper <- function(alpha, Z, size, b, q){
  CIobj = function(theta, alpha, Z, size, b, q){
    return((pvalRight(Z = Z, size = size, theta = theta, b = b, q = q) - alpha)^2)
  }
  U = stats::optim(par = .5, fn = CIobj,
                   alpha = 1-alpha, Z = Z, size = size, b = b, q = q,
                   method = "Brent", lower = 0, upper = 1)
  return(U$par)
}

