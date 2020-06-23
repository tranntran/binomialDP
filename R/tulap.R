#' @title The Truncated-Uniform-Laplace (Tulap) Distribution
#' @name tulap
#' @aliases ptulap
#' @aliases rtulap
#'
#' @description Distribution function and random generation for the Tulap
#'   distribution.
#' @param t vector of quantiles
#' @param m vector of medians
#' @param b vector of Discrete Laplace noise parameters, obtained by
#'   \eqn{exp(-\epsilon)}
#' @param q vector of truncated quantiles
#' @param n number of observations
#'
#' @examples
#' sample <- rtulap(1000, 1, 0.3, 0.05)
#' sample[1:10]
#' plot(density(sample))
#' ptulap(0.5, 0, 0.3, 0.05)
#'
#' @return \code{ptulap} gives the distribution function and \code{rtulap}
#'   generates random derviates.
#' @references \code{ptulap} and \code{rtulap} are based on Awan, Jordan
#'   Alexander, and Aleksandra Slavkovic. 2020. "Differentially Private
#'   Inference for Binomial Data". Journal of Privacy and Confidentiality 10
#'   (1). \url{https://doi.org/10.29012/jpc.725}.
NULL

#' @rdname tulap
#' @export
#'
rtulap <- function (n, m = 0, b = 0, q = 0) {
  if(q >= 0){
    alpha = .95
    lcut = q/2
    rcut = q/2

    # Approximates M such that given n Bernoulli trials with success rate prob,
    # is such that alpha of the times, there are at least n successes among M. n
    # - the number of trials that we want more than success of prob - Bernoulli
    # probability success of each iid trial alpha - probability that among the M
    # trials, at least n successes.
    approx.trials <- function (n, prob = 1, alpha = 0) {
      # Solve a quadratic form for this:
      a = prob^2
      b = -((2 * n * prob) + ((stats::qnorm(alpha)^2) * prob * (1 - prob)))
      c = n^2
      return ((-b + sqrt(b^2 - (4 * a * c))) / (2 * a))
    }

    # Calculate actual amount needed
    q = lcut + rcut
    n2 = approx.trials(n, prob=(1 - q), alpha=alpha)

    # Sample from the original Tulambda distribution
    geos1 = stats::rgeom(n2, prob=(1 - b))
    geos2 = stats::rgeom(n2, prob=(1 - b))
    unifs = stats::runif(n2, min=(-1/2), max=(1/2))
    samples = m + geos1 - geos2 + unifs

    # Cut the tails based on the untampered CDF (ie no cuts)
    probs = ptulap(samples, m = m, b = b)
    is.mid = (lcut <= probs) & (probs <= (1 - rcut))

    # Abuse the NA property of R wrt arithmetics
    mids = samples[is.mid]
    while ({len = base::length(mids); len} < n) {
      diff = n - len
      mids = c(mids, rtulap(diff, m=m, b=b,q=q))
    }
    return (mids[1:n])
  }
  geos1 = stats::rgeom(n2, prob=(1 - b))
  geos2 = stats::rgeom(n2, prob=(1 - b))
  unifs = stats::runif(n2, min=(-1/2), max=(1/2))
  samples = m + geos1 - geos2 + unifs
  return(samples)
}


#' @rdname tulap
#' @export
ptulap <- function (t, m = 0, b = 0, q = 0) {
  lcut = q/2
  rcut = q/2
  # Normalize
  t = t - m

  # Split the positive and negsative t calculations, and factor out stuff
  r = base::round(t)
  g = -base::log(b)
  l = base::log(1 + b)
  k = 1 - b
  negs = base::exp((r * g) - l + base::log(b + ((t - r + (1/2)) * k)))
  poss = 1 - base::exp((r * (-g)) - l + base::log(b + ((r - t + (1/2)) * k)))

  # Check for infinities
  negs[base::is.infinite(negs)] <- 0
  poss[base::is.infinite(poss)] <- 0

  # Truncate wrt the indicator on t's positivity
  is.leq0 = t <= 0
  trunc = (is.leq0 * negs) + ((1 - is.leq0) * poss)

  # Handle the cut adjustment and scaling
  q = lcut + rcut
  is.mid = (lcut <= trunc) & (trunc <= (1 - rcut))
  is.rhs = (1 - rcut) < trunc
  return (((trunc - lcut) / (1 - q)) * is.mid + is.rhs)
}
