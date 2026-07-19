#' @title Simulate piecewise exponential time-to-event outcomes
#'
#' @description Simulate time-to-event outcomes using the piecewise constant
#'   hazard exponential function.
#'
#' @param n integer. The number of random samples to generate. Default is
#'   `n = 1`.
#' @param hazard vector. Finite non-negative constant hazard rates for
#'   exponential failures. If the final rate is zero, `maxtime` must be
#'   supplied so that subjects without an event can be administratively censored.
#' @param cutpoints finite, positive, strictly increasing vector of interior
#'   times at which the hazard rate changes. The number of hazard rates must be
#'   one greater than the number of cutpoints. Use `NULL` for a constant hazard.
#' @param maxtime scalar. Optional administrative censoring time. When
#'   supplied, it must be later than every cutpoint.
#'
#' @details See [pwe_impute()] for details.
#'
#' @return A data frame with simulated follow-up times (`time`) and respective
#'   event indicator (`event`, 1 = event occurred, 0 = censoring).
#'
#' @importFrom PWEALL rpwe qpwe pwe
#' @importFrom stats rexp
#' @importFrom utils tail
#' @export
#'
#' @examples
#' pwe_sim(10, hazard = c(0.005, 0.001), cutpoints = 3, maxtime = 36)
#' y <- pwe_sim(n = 1, hazard = c(2.585924e-02, 3.685254e-09),
#'              cutpoints = 12)
pwe_sim <- function(n = 1, hazard = 1, cutpoints = NULL, maxtime = NULL) {
  validate_cutpoints(cutpoints)
  validate_piecewise_hazard(hazard, cutpoints)
  validate_maxtime(maxtime, cutpoints)

  if (is.null(maxtime) && tail(hazard, 1) == 0) {
    stop("'maxtime' must be supplied when the final hazard rate is zero")
  }

  if (length(hazard) == 1) {
    # stats::rexp(rate = 0) returns NaN, whereas a zero constant hazard
    # represents no event before administrative censoring.
    ret <- if (hazard == 0) rep(Inf, n) else rexp(n, rate = hazard)
  } else {
    ret <- PWEALL::rpwe(n, rate = hazard, tchange = c(0, cutpoints))$r
  }

  if (!is.null(maxtime)) {
    min_time <- pmin(ret, maxtime)
    event <- as.numeric(ret == min_time)
    dat <- data.frame(time = min_time, event = event)
  } else {
    dat <- data.frame(time = ret, event = rep(1, length(ret)))
  }

  return(dat)
}


#' @title Impute piecewise exponential time-to-event outcomes
#'
#' @description Imputation of time-to-event outcomes using the piecewise
#'   constant hazard exponential function conditional on observed exposure.
#'
#' @inheritParams pwe_sim
#' @param time vector. The observed time for patient that have had no event or
#'   passed `maxtime`.
#'
#' @details If a subject is event-free at time \eqn{s < t}, then the conditional
#'   probability is
#'
#'   \deqn{F_{T | s}(t | s) = P(T \le t | T > s) = \frac{F(t) - F(s)}{1 - F(s)}}
#'
#'   where \eqn{F(\cdot)} is the cumulative distribution function of the
#'   piecewise exponential (PWE) distribution. Equivalently, \eqn{F(t) = 1 -
#'   S(t)}, where `S(t)` is the survival function. If \eqn{U \sim Unif(0, 1)},
#'   then we can generate an event time (conditional on being event free up
#'   until \eqn{s}) as
#'
#'   \deqn{F^{-1}(U(1 - F(s)) + F(s))}
#'
#'   If \eqn{s = 0}, this is equivalent to a direct unconditional sample from
#'   the PWE distribution.
#'
#' @return A data frame with simulated follow-up times (`time`) and respective
#'   event indicator (`event`, 1 = event occurred, 0 = censoring).
#' @export
#'
#' @examples
#' pwe_impute(time = c(3, 4, 5), hazard = c(0.002, 0.01), cutpoints = 12)
#' pwe_impute(time = c(3, 4, 5), hazard = c(0.002, 0.01), cutpoints = 12,
#'            maxtime = 36)
#' pwe_impute(time = 19.621870008, hazard = c(2.585924e-02, 3.685254e-09),
#'            cutpoints = 12, maxtime = 36)
pwe_impute <- function(time, hazard, cutpoints = NULL, maxtime = NULL) {
  validate_cutpoints(cutpoints)
  validate_piecewise_hazard(hazard, cutpoints)
  validate_maxtime(maxtime, cutpoints)

  # Check: 'time' is positive integer
  if (any(time < 0)) {
    stop("'time' must always be positive")
  }

  if (!is.null(maxtime)) {
    if (any(maxtime < time)) {
      stop(
        "'maxtime' must be greater than or equal to all observed 'time' values"
      )
    }
  }

  if (is.null(maxtime) && tail(hazard, 1) == 0) {
    stop("'maxtime' must be supplied when the final hazard rate is zero")
  }

  # Use inverse CDF to get conditional samples
  interval_starts <- c(0, cutpoints)
  Fs <- PWEALL::pwe(t = time, rate = hazard, tchange = interval_starts)$dist
  U <- runif(length(time))
  time_imp <- PWEALL::qpwe(U * (1 - Fs) + Fs, hazard, interval_starts)$q

  # Check: impute timed occur after landmark observed times
  if (any(time > time_imp)) {
    stop("Imputed times cannot precede the observed times")
  }

  if (!is.null(maxtime)) {
    min_time <- pmin(time_imp, maxtime)
    event <- as.numeric(time_imp == min_time)
    dat <- data.frame(time = min_time, event = event)
  } else {
    dat <- data.frame(time = time_imp, event = rep(1, length(time)))
  }

  return(dat)
}


#' @title Cumulative distribution function of the PWE for a vectorized hazard
#'   rate parameter
#'
#' @description Extends [PWEALL::pwe()] to allow for vectorization over the
#'   hazard rates.
#'
#' @param hazard matrix. A matrix of hazard rate parameters with number of
#'   columns one greater than the length of the `cutpoints` vector. The number
#'   of rows can be anything, and is typically dictated by the number of MCMC
#'   draws.
#' @param end_of_study finite positive time at which the cumulative event
#'   probability is evaluated. It must be greater than every cutpoint.
#' @inheritParams pwe_sim
#' @inheritParams survival_adapt
#'
#' @return A vector of (0, 1) probabilities from evaluation of the PWE
#'   cumulative distribution function. Length of the vector matches the number
#'   of rows of the `hazard` matrix parameter.
#'
#' @export
ppwe <- function(hazard, end_of_study, cutpoints = NULL) {
  validate_cutpoints(cutpoints)
  validate_endpoint_time(end_of_study, cutpoints, "end_of_study")
  validate_hazard_matrix(hazard, cutpoints)

  interval_lower <- c(0, cutpoints)
  interval_upper <- c(cutpoints, Inf)
  duration <- pmax(0, pmin(end_of_study, interval_upper) - interval_lower)
  1 - exp(-drop(hazard %*% duration))
}
