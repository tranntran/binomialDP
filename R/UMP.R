#' @title Calculating Simple and One-Sided DP-UMP Tests
#' @name UMP
#' @aliases ump
#' @aliases umpleft
#' @aliases umpright
#'
#' @param theta The success probability for each trial (parameter \eqn{\theta}
#'   in Binomial(n, \eqn{\theta}))
#' @param size The number of trials in Binomial distribution (parameter n in
#'   Binomial(n, \eqn{\theta}))
#' @param alpha Level of the tests
#' @param epsilon Parameter \eqn{\epsilon} in \eqn{(\epsilon, \delta)}-DP
#' @param delta Parameter \eqn{\delta} in \eqn{(\epsilon, \delta)}-DP
#'
#' @return A vector of one-sided DP-UMP tests
#' @examples
#' #leftUMP
#' left <- umpLeft(theta = 0.4, size = 10, alpha = 0.05, epsilon = 1, delta = 0.01)
#'
#' #rightUMP
#' right <- umpRight(theta = 0.4, size = 10, alpha = 0.05, epsilon = 1, delta = 0.01)
#'
#' #plot
#' plot(left, type = "l", main = "One-sided DP UMP",
#'      xlab = "x", ylab = "Phi (Probability of rejecting the null hypothesis)")
#' lines(right, lty = 2, col = "blue")
#' legend("left", legend=c("Left UMP", "Right UMP"), col=c("black", "blue"), lty=1:2, lwd = 1.5)
#'
#' @references Awan, Jordan Alexander, and Aleksandra Slavkovic. 2020.
#'   "Differentially Private Inference for Binomial Data". Journal of Privacy
#'   and Confidentiality 10 (1). \url{https://doi.org/10.29012/jpc.725}.
#' @seealso Calculating unbiased two-sided DP-UMP tests (\code{\link{UMPU}}) and
#'   asymptotically unbiased two-sided DP-UMP tests (\code{\link{umpuApprox}})
NULL

#' @rdname UMP
#' @export
#'
umpLeft <- function(theta, size, alpha, epsilon, delta){
  b = base::exp(-epsilon)
  q = 2*delta*b/(1-b+2*delta*b)
  values = base::seq(0, size)
  B = stats::dbinom(values, size = size, prob = theta)

  obj = function(s){
    phi = ptulap(t = values-s, m = 0, b = b, q = q)
    return(B%*%phi - alpha)
  }

  lower = -1
  while(obj(lower) < 0)
    lower = lower*2
  upper = 1
  while(obj(upper) > 0)
    upper = upper*2
  result = stats::uniroot(f = obj, interval = c(lower, upper))
  s = result$root
  phi = ptulap(t = values-s, m = 0, b = b, q = q)
  return(phi)
}

#' @rdname UMP
#' @export

umpRight <- function(theta, size, alpha, epsilon, delta){
  b = base::exp(-epsilon)
  q = 2*delta*b/(1-b+2*delta*b)
  values = base::seq(0, size)
  B = stats::dbinom(values, size = size, prob = theta)

  obj = function(s){
    phi = ptulap(t = values-s, m = 0, b = b, q = q)
    return(B%*%phi - alpha)
  }

  lower = -1
  while(obj(lower) < 0)
    lower = lower*2
  upper = 1
  while(obj(upper) > 0)
    upper = upper*2
  result = stats::uniroot(f = obj, interval = c(lower, upper))
  s = result$root
  phi = ptulap(t = values-s, m = 0, b = b, q = q)
  return(1-phi)
}
