#' @title Derive piecewise-constant hazard rates from cumulative event
#'   probabilities
#'
#' @description Converts cumulative event probabilities at specified follow-up
#'   times into the corresponding piecewise-constant hazard rates. This is
#'   useful for expressing data-generation assumptions in terms of clinically
#'   interpretable event probabilities.
#'
#' @param probs A required numeric vector of finite cumulative event
#'   probabilities in `[0, 1)` at
#'   each cutpoint and at `endtime`, in that order. Its length must be one
#'   greater than the number of cutpoints. With no cutpoints, supply a single
#'   probability at `endtime`. Values must be non-decreasing and are not
#'   recycled.
#' @param cutpoints `NULL` (the default), or a numeric vector of finite,
#'   positive, strictly increasing interior
#'   times at which the event hazard changes. `NULL`
#'   corresponds to a simple (non-piecewise) exponential model.
#' @param endtime A required single finite, positive numeric value giving the
#'   follow-up time corresponding to the final element of `probs`. It must be
#'   later than every cutpoint and use the same time unit.
#'
#' @details Given \eqn{J-1} interior cutpoints, then there are J intervals
#'   defined as: \eqn{[s_0, s_1)}, \eqn{[s_1, s_2)}, \eqn{\dots}, \eqn{[s_{J-1},
#'   s_{J})}, with conditions that \eqn{s_0 = 0} and \eqn{s_J = \infty}. Each
#'   interval corresponds to constant hazard \eqn{\lambda_j}. This is the
#'   PWEALL representation of the continuous generating hazard. Changing the
#'   value at an isolated cutpoint does not alter the cumulative probabilities
#'   calculated here. When observed event times are assigned to analysis
#'   intervals, `goldilocks` uses `(s_{j-1}, s_j]`, matching the survival
#'   counting-process convention, so an event at \eqn{s_j} belongs to the
#'   interval ending there.
#'
#' @return A numeric vector of non-negative hazard rates, with one value for
#'   each interval defined by `cutpoints` and `endtime`.
#'
#' @export
#'
#' @examples
#' lambda <- prop_to_haz(0.15, endtime = 36) # 15% probability at 36-months
#' all.equal(pexp(36, lambda), 0.15)
#'
#' # 15% probability at 12-months, and 30% at 24-months
#' prop_to_haz(c(0.15, 0.30), 12, 24)
#' PWEALL::pwe(12, prop_to_haz(c(0.15, 0.30), 12, 24), c(0, 12))$dist
#' PWEALL::pwe(24, prop_to_haz(c(0.15, 0.30), 12, 24), c(0, 12))$dist
prop_to_haz <- function(probs, cutpoints = NULL, endtime) {
  validate_probability_vector(probs, "probs", upper_open = TRUE)

  if (any(diff(probs) < 0)) {
    stop("'probs' must be non-decreasing over time")
  }

  validate_cutpoints(cutpoints)
  validate_endpoint_time(endtime, cutpoints, "endtime")

  J <- length(cutpoints) + 1L
  if (length(probs) != J) {
    stop("Length of 'probs' must be one greater than length of 'cutpoints'")
  }

  lambda <- vector(length = J)
  cumulative_hazard <- probability_to_cumulative_hazard(probs)

  if (J == 1) {
    lambda <- cumulative_hazard / endtime
  } else {
    s <- c(0, cutpoints, endtime)
    s_diff <- diff(s)
    lambda[1] <- cumulative_hazard[1] / s_diff[1]
    for (j in 2:J) {
      offset <- sum(lambda[1:(j - 1)] * s_diff[1:(j - 1)])
      lambda[j] <- (cumulative_hazard[j] - offset) / s_diff[j]
    }
  }

  return(lambda)
}
