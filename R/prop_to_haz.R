#' @title Estimate plausible piecewise constant hazard rates from summary
#'   summary event proportions
#'
#' @description Given estimates of the event probability at one or more fixed
#'   times, the corresponding piecewise hazard rates can be determined through
#'   closed-form formulae. This utility function can be useful when simulating
#'   trial datasets with plausible event rates.
#'
#' @param probs vector. Probabilities of the event (i.e. cumulative incidence
#'   probabilities) at one or more time point. If only a single value is given,
#'   then it is assumed that this is the probability at the \code{endtime}.
#' @param cutpoints vector. Times at which the baseline hazard changes. Default
#'   is \code{cutpoints = 0}, which corresponds to a simple (non-piecewise)
#'   exponential model.
#' @param endtime scalar. Time at which final element in \code{probs}
#'   corresponds to. Typically this would be the study endpoint time.
#'
#' @details Given \eqn{J-1} internal cut-points, then there are J intervals
#'   defined as: \eqn{[s_0, s_1)}, \eqn{[s_1, s_2)}, \eqn{\dots}, \eqn{[s_{J-1},
#'   s_{J})}, with conditions that \eqn{s_0 = 0} and \eqn{s_J = \infty}. Each
#'   interval corresponds to constant hazard \eqn{\lambda_j}.
#'
#' @return Vector of constant hazard rates for each time piece defined by
#'   `cutpoints`.
#'
#' @export
#'
#' @examples
#' lambda <- prop_to_haz(0.15, endtime = 36) # 15% probability at 36-months
#' all.equal(pexp(36, lambda), 0.15)
#'
#' # 15% probability at 12-months, and 30% at 24-months
#' prop_to_haz(c(0.15, 0.30), c(0, 12), 24)
#' PWEALL::pwe(12, prop_to_haz(c(0.15, 0.30), c(0, 12), 24), c(0, 12))$dist
#' PWEALL::pwe(24, prop_to_haz(c(0.15, 0.30), c(0, 12), 24), c(0, 12))$dist
prop_to_haz <- function(probs, cutpoints = 0, endtime) {

  if (cutpoints[1] != 0) {
    stop("First element of 'cutpoints' should be 0")
  }

  J <- length(cutpoints)
  lambda <- vector(length = J)

  if (J == 1) {
    lambda <- -log(1 - probs) / endtime
  } else {
    s <- c(cutpoints, endtime)
    s_diff <- diff(s)
    lambda[1] <- -log(1 - probs[1]) / s_diff[1]
    for (j in 2:J) {
      offset <- sum(lambda[1:(j-1)] * s_diff[1:(j-1)])
      lambda[j] <- (-log(1 - probs[j]) - offset) / s_diff[j]
    }
  }

  return(lambda)

}
