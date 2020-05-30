#' Asymptotically Unbiased Two-Sided DP-UMP Tests
#' @aliases umpuapprox
#'
#' @param theta The success probability for each trial (parameter \eqn{\theta}
#'   in Binomial(n, \eqn{\theta}))
#' @param size The number of trials in Binomial distribution (parameter n in
#'   Binomial(n, \eqn{\theta}))
#' @param alpha Level of the tests
#' @param epsilon Parameter \eqn{\epsilon} in \eqn{(\epsilon, \delta)}-DP
#' @param delta Parameter \eqn{\delta} in \eqn{(\epsilon, \delta)}-DP
#'
#' @return A vector of asymptotically unbiased two-sided DP-UMP tests
#' @export
#' @references Awan, Jordan Alexander, and Aleksandra Slavkovic. 2020.
#'   "Differentially Private Inference for Binomial Data". Journal of Privacy
#'   and Confidentiality 10 (1). \url{https://doi.org/10.29012/jpc.725}.
#' @seealso Calculating simple and one-sided DP-UMP tests (\code{\link{umpLeft}}
#'   or \code{\link{umpRight}}) and unbiased two-sided DP-UMP tests
#'   (\code{\link{UMPU}})
#' @examples
#' #Comparing unbiased DP-UMP tests obtained by umpuApprox and UMPU
#' asymp <- umpuApprox(theta = 0.4, size = 10, alpha = 0.05, epsilon = 1, delta = 0.01)
#' twoside <- UMPU(theta = 0.4, size = 10, alpha = 0.05, epsilon = 1, delta = 0.01)
#'
#' #Plot the probability of rejecting the null hypothesis based on x
#' plot(asymp, type = "l", lwd = 1.5, col = "red", xlab = "x",
#'   ylab = "Phi (Probability of rejecting the null hypothesis)",
#'   main = "Unbiased DP-UMP Tests")
#' lines(twoside, type = "l", lwd = 1.5, lty = 2, col = "blue")
#' legend("topleft", legend=c("Asymptotically UMPU", "UMPU"),
#'   col=c("red", "blue"), lty=1:2, lwd = 1.5)
#'
umpuApprox <-  function(theta, size, alpha, epsilon, delta){
  b = base::exp(-epsilon)
  q = 2*delta*b/(1-b+2*delta*b)
  values = base::seq(0, size)
  B = stats::dbinom(values, size = size, prob = theta)
  BX = B*(values-size*theta)
  k = size*theta
  greaterK = (values >= k)

  obj = function(s){
    F1 = ptulap(t = values- k-s, m = 0, b = b, q = q)
    F2 = ptulap(t = k-values-s, m = 0, b = b, q = q)
    phi = F1*greaterK + F2*(1-greaterK)
    return(B%*%phi-alpha)
  }
  lower = -1
  while(obj(lower) < 0)
    lower = lower*2
  upper = 1
  while(obj(upper) > 0)
    upper = upper*2

  result = stats::uniroot(f = obj,interval = c(lower, upper))
  s = result$root

  F1 = ptulap(t = values-k-s, m = 0, b = b, q = q)
  F2 = ptulap(t = k-values-s, m = 0, b = b, q = q)
  phi = F1*greaterK + F2*(1-greaterK)
  return(phi)
}
