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
#'   is \code{NULL}, which corresponds to a simple (non-piecewise) exponential
#'   model.
#' @param endtime scalar. Time at which final element in \code{probs}
#'   corresponds to. Typically this would be the study endpoint time.
#'
#' @details Given \eqn{J-1} cut-points, then there are J intervals defined as:
#'   \eqn{[s_0, s_1)}, \eqn{[s_1, s_2)}, \eqn{\dots}, \eqn{[s_{J-1}, s_{J})},
#'   with conditions that \eqn{s_0 = 0} and \eqn{s_J = \infty}. Each interval
#'   corresponds to constant hazard \eqn{\lambda_j}.
#'
#' @return Vector of constant hazard rates.
#'
#' @export
#'
#' @examples
#' lambda <- haz_est(0.15, endtime = 36) # 15% probability at 36-months
#' all.equal(pexp(36, lambda), 0.15)
#'
#'
#' haz_est(c(0.15, 0.30), 12, 24) # 15% probability at 12-months, and 30% at
#' PWEALL::pwe(12, haz_est(c(0.15, 0.30), 12, 24), c(0, 12))$dist
#' PWEALL::pwe(24, haz_est(c(0.15, 0.30), 12, 24), c(0, 12))$dist
haz_est <- function(probs, cutpoints = NULL, endtime) {

  if (!is.null(cutpoints)) {
    if (cutpoints[1] == 0) {
      stop("First element of 'cutpoints' should not be 0")
    }
    J <- length(cutpoints) + 1
  } else {
    J <- 1
  }

  lambda <- vector(length = J)

  if (J == 1) {
    lambda <- -log(1 - probs) / endtime
  } else {
    s <- c(0, cutpoints, endtime)
    s_diff <- diff(s)
    lambda[1] <- -log(1 - probs[1]) / s_diff[1]
    for (j in 2:J) {
      offset <- sum(lambda[1:(j-1)] * s_diff[1:(j-1)])
      lambda[j] <- (-log(1 - probs[j]) - offset) / s_diff[j]
    }
  }

  return(lambda)

}
