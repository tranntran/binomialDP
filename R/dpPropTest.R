#' DP Population Proportion Test

#' @param Z A Binomial sample with Tulap noise
#' @param size The number of trials in Binomial distribution (parameter n in
#'   Binomial(n, \eqn{\theta}))
#' @param theta The success probability for each trial (parameter \eqn{\theta}
#'   in Binomial(n, \eqn{\theta}))
#' @param b Discrete Laplace noise parameters, obtained by \eqn{exp(-\epsilon)}
#' @param q The truncated quantiles
#' @param alternative The type of alternative hypothesis, must be
#'   \code{"two.sided"}, \code{"less"}, or \code{"greater"}
#' @param alpha The confidence level of the confidence interval
#'
#' @return A table summary of the sample estimate, the p-value of the test, and
#'   the confidence interval of \eqn{\theta}
#' @export
#'
#' @examples dpPropTest(9, 30, 0.5, exp(-1), 0.05, "less", 0.05)
dpPropTest <- function(Z, size, theta, b, q, alternative, alpha){
  if(alternative == "two.sided"){
    pval = pvalTwoSide(Z, size, theta, b, q)
    lowerLimit = CITwoSide(alpha, Z, size, b, q)[1]
    upperLimit = CITwoSide(alpha, Z, size, b, q)[2]
    h = paste("Ha: True theta not equal to", theta)
  } else if(alternative == "less"){
    pval = pvalLeft(Z, size, theta, b, q)
    lowerLimit = 0
    upperLimit = CIUpper(alpha, Z, size, b, q)
    h = paste("Ha: True theta less than", theta)
  } else if(alternative == "greater"){
    pval = pvalRight(Z, size, theta, b, q)
    lowerLimit = CILower(alpha, Z, size, b, q)
    upperLimit = 1
    h = paste("Ha: True theta greater than", theta)
  }
  result = rbind(Z, size, Z/size, pval, lowerLimit, upperLimit)
  rownames(result) = c("Data", "Size", "Sample Estimates", "p-value", "Lower Limit", "Upper Limit")
  colnames(result) = h
  return(result)
}
