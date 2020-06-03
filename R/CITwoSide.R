#' Finding Asymptotically Unbiased Two-Sided Confidence Intervals
#' @name CITwoSide
#' @aliases citwoside
#'
#' @param alpha The confidence level of the returned confidence interval
#' @param Z A Binomial sample with Tulap noise
#' @param size The number of trials in Binomial distribution (parameter n in
#'   Binomial(n, \eqn{\theta}))
#' @param b The discrete Laplace noise parameter, obtained by
#'   \eqn{exp(-\epsilon)}
#' @param q The truncated quantiles
#'
#' @return An asymptotically unbiased two-sided confindence interval
#' @export
#' @references Awan, Jordan Alexander, and Aleksandra Slavkovic. 2020.
#'   "Differentially Private Inference for Binomial Data". Journal of Privacy
#'   and Confidentiality 10 (1). \url{https://doi.org/10.29012/jpc.725}.
#' @seealso Finding one-sided confidence intervals \code{\link{CILower}} and
#'   \code{\link{CIUpper}}
#'
#' @examples
#' CITwoSide(alpha = 0.05, Z = 18, size = 30, b = exp(-1), q = 0.05) #(0.411, 0.767)
#'
CITwoSide <- function(alpha , Z, size, b, q){
  mle = Z/size
  mle = base::max(base::min(mle, 1), 0)
  CIobj2 = function(theta, alpha, Z, size, b, q){
    return((pvalTwoSide(Z = Z, theta = theta, size = size, b = b, q = q) - alpha)^2)
  }

  if(mle > 0){
    L = stats::optim(par = mle/2, fn = CIobj2,
                     alpha = alpha, Z = Z, size = size, b = b, q = q,
                     method = "Brent", lower = 0, upper = mle)
    L = L$par
  } else {
    L = 0
  }

  if(mle < 1){
    U = stats::optim(par = (1-mle)/2, fn = CIobj2,
                     alpha = alpha, Z = Z, size = size, b = b, q = q,
                     method = "Brent", lower = mle, upper = 1)
    U = U$par
  } else {
    U = 1
  }

  CI = c(L, U)
  return(CI)
}
