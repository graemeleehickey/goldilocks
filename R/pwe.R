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
#' @param cutpoints vector. The change-point vector indicating time when the
#'   hazard rates change. Note the first element of `cutpoints` should
#'   always be 0.
#' @param maxtime scalar. Maximum time before end of study.
#'
#' @details See [pwe_impute()] for details.
#'
#' @return A data frame with simulated follow-up times (`time`) and respective
#'   event indicator (`event`, 1 = event occurred, 0 = censoring).
#'
#' @importFrom PWEALL rpwe qpwe pwe
#' @importFrom stats rexp
#' @export
#'
#' @examples
#' pwe_sim(10, hazard = c(0.005, 0.001), cutpoints = c(0, 3), maxtime = 36)
#' y <- pwe_sim(n = 1, hazard = c(2.585924e-02, 3.685254e-09),
#'              cutpoints = c(0, 12))
pwe_sim <- function(n = 1, hazard = 1, cutpoints = 0, maxtime = NULL) {
  validate_piecewise_hazard(hazard, cutpoints)

  # Check: first element of 'cutpoints' should be 0
  if (cutpoints[1] != 0) {
    stop("First element of 'cutpoints' should be 0")
  }

  # Check: 'cutpoints' is increasing
  if (is.unsorted(cutpoints)) {
    stop("'cutpoints' should be in increasing order")
  }

  # Check: 'maxtime' is positive or NULL
  if (!is.null(maxtime)) {
    if (maxtime <= 0 | length(maxtime) > 1) {
      stop("'maxtime' must be a postive single value")
    }
  }

  if (is.null(maxtime) && tail(hazard, 1) == 0) {
    stop("'maxtime' must be supplied when the final hazard rate is zero")
  }

  if (length(hazard) == 1) {
    # stats::rexp(rate = 0) returns NaN, whereas a zero constant hazard
    # represents no event before administrative censoring.
    ret <- if (hazard == 0) rep(Inf, n) else rexp(n, rate = hazard)
  } else {
    ret <- PWEALL::rpwe(n, rate = hazard, tchange = cutpoints)$r
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
#' pwe_impute(time = c(3, 4, 5), hazard = c(0.002, 0.01), cutpoints = c(0, 12))
#' pwe_impute(time = c(3, 4, 5), hazard = c(0.002, 0.01), cutpoints = c(0, 12),
#'            maxtime = 36)
#' pwe_impute(time = 19.621870008, hazard = c(2.585924e-02, 3.685254e-09),
#'            cutpoints = c(0, 12), maxtime = 36)
pwe_impute <- function(time, hazard, cutpoints = 0, maxtime = NULL) {
  validate_piecewise_hazard(hazard, cutpoints)

  # Check: 'time' is positive integer
  if (any(time < 0)) {
    stop("'time' must always be positive")
  }

  # Check: 'maxtime' is positive or NULL
  if (!is.null(maxtime)) {
    if (maxtime <= 0 | length(maxtime) > 1) {
      stop("'maxtime' must be a postive single value")
    }
    if (any(maxtime < time)) {
      stop("'maxtime' must be greater than or equal to all observed 'time' values")
    }
  }

  if (is.null(maxtime) && tail(hazard, 1) == 0) {
    stop("'maxtime' must be supplied when the final hazard rate is zero")
  }

  # Use inverse CDF to get conditional samples
  Fs <- PWEALL::pwe(t = time, rate = hazard, tchange = cutpoints)$dist
  U <- runif(length(time))
  time_imp <- PWEALL::qpwe(U * (1 - Fs) + Fs, hazard, cutpoints)$q

  # impute1 <- function(s) {
  #   Fs <- bayesDP::ppexp(s, hazard, cutpoints)
  #   u <- runif(1)
  #   PWEALL::qpwe(u*(1 - Fs) + Fs, hazard, cutpoints)$q
  #   #msm::qpexp(u*(1 - Fs) + Fs, hazard, cutpoints)
  #   #rpact::qpwexp(u*(1 - Fs) + Fs, lambda = hazard, s = cutpoints)
  # }
  #
  # time_imp <- sapply(time, impute1)

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
#'   columns equal to the length of the `cutpoints` vector. The number of
#'   rows can be anything, and is typically dictated by the number of MCMC
#'   draws.
#' @inheritParams pwe_sim
#' @inheritParams survival_adapt
#'
#' @return A vector of (0, 1) probabilities from evaluation of the PWE
#'   cumulative distribution function. Length of the vector matches the number
#'   of rows of the `hazard` matrix parameter.
#'
#' @export
ppwe <- function(hazard, end_of_study, cutpoints) {
  validate_hazard_matrix(hazard, cutpoints)

  interval_upper <- c(cutpoints[-1], Inf)
  duration <- pmax(0, pmin(end_of_study, interval_upper) - cutpoints)
  1 - exp(-drop(hazard %*% duration))
}
