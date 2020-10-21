#' @title Simulate piecewise exponential time-to-event outcomes
#'
#' @description Simulate time-to-event outcomes using the piecewise constant
#'   hazard exponential function.
#'
#' @param n integer. The number of random samples to generate. Default is
#'   \code{n=1}.
#' @param hazard vector. The constant hazard rates for exponential failures.
#' @param cutpoint vector. The change-point vector indicating time when the
#'   hazard rates change. Note the first element of \code{cutpoint} should
#'   always be 0.
#' @param maxtime scalar. Maximum time before end of study.
#'
#' @details See \code{\link{pwe_impute}} for details.
#'
#' @return A data frame with simulated follow-up times (\code{time}) and
#'   respective event indicator (\code{event}, 1 = event occurred, 0 =
#'   censoring).
#'
#' @importFrom PWEALL rpwe qpwe pwe
#' @export
#'
#' @examples
#' pwe_sim(10, hazard = c(0.005, 0.001), cutpoint = c(0, 3), maxtime = 36)
#' y <- pwe_sim(n = 1, hazard = c(2.585924e-02, 3.685254e-09),
#'              cutpoint = c(0, 12))
pwe_sim <- function(n = 1, hazard = 1, cutpoint = 0, maxtime = NULL) {

  # Check: 'cutpoint' should be same length as hazard
  if (length(cutpoint) != length(hazard)) {
    stop("length of t must be equal to length of rate")
  }

  # Check: first element of 'cutpoint' should be 0
  if (cutpoint[1] != 0) {
    stop("first element of t should be 0")
  }

  # Check: 'cutpoint' is increasing
  if (is.unsorted(cutpoint)) {
    stop("t should be in increasing order")
  }

  # Check: 'maxtime' is positive or NULL
  if (!is.null(maxtime)) {
    if (maxtime <= 0 | length(maxtime) > 1) {
      stop("The 'maxtime' must be a postive single value!")
    }
  }

  ret <- PWEALL::rpwe(n, rate = hazard, tchange = cutpoint)$r

  if (!is.null(maxtime)) {
    min_time <- pmin(ret, maxtime)
    event <- as.numeric(ret == min_time)
    dat <- data.frame(time = min_time, event = event)
  }
  else {
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
#'   passed \code{maxtime}.
#'
#' @details If a subject is event-free at time \eqn{s < t}, then the conditional
#'   probability \eqn{F_{T \| s}|(t \| s) = P[T \le \| T > s] = (F(t) - F(s)) /
#'   (1 - F(s))}, where \eqn{F(\cdot)} is the cumulative distribution function
#'   of the piecewise exponential (PWE) distribution. Equivalently, \eqn{F(t) =
#'   1 - S(t)}, where \code{S(t)} is the survival function. If \eqn{U \sim
#'   Unif(0, 1)}, then we can generate an event time (conditional on being event
#'   free up until \eqn{s}) as \eqn{F^{-1}(U(1-F(s)) + F(s))}. Note: if \eqn{s =
#'   0}, then this is the equivalent of a direct (unconditional) sample from the
#'   PWE distribution.
#'
#' @return A data frame with simulated follow-up times (\code{time}) and
#'   respective event indicator (\code{event}, 1 = event occurred, 0 =
#'   censoring).
#' @importFrom bayesDP ppexp
#' @export
#'
#' @examples
#' pwe_impute(time = c(3, 4, 5), hazard = c(0.002, 0.01), cutpoint = c(0, 12))
#' pwe_impute(time = c(3, 4, 5), hazard = c(0.002, 0.01), cutpoint = c(0, 12),
#'            maxtime = 36)
#' pwe_impute(time = 19.621870008, hazard = c(2.585924e-02, 3.685254e-09),
#'            cutpoint = c(0, 12), maxtime = 36)
pwe_impute <- function(time, hazard, cutpoint = 0, maxtime = NULL) {

  # Check: 'hazard' is positive
  if (any(hazard < 0)) {
    stop("At least one of the hazard rate(s) is less than 0!")
  }

  # Check: 'time' is positive integer
  if (any(time < 0)) {
    stop("The time has to always be positive!")
  }

  # Check: 'maxtime' is positive or NULL
  if (!is.null(maxtime)) {
    if (maxtime <= 0 | length(maxtime) > 1) {
      stop("The 'maxtime' must be a postive single value!")
    }
  }

  # Use inverse CDF to get conditional samples
  Fs <- PWEALL::pwe(t = time, rate = hazard, tchange = cutpoint)$dist
  U <- runif(length(time))
  time_imp <- PWEALL::qpwe(U*(1 - Fs) + Fs, hazard, cutpoint)$q

  # impute1 <- function(s) {
  #   Fs <- bayesDP::ppexp(s, hazard, cutpoint)
  #   u <- runif(1)
  #   PWEALL::qpwe(u*(1 - Fs) + Fs, hazard, cutpoint)$q
  #   #msm::qpexp(u*(1 - Fs) + Fs, hazard, cutpoint)
  #   #rpact::qpwexp(u*(1 - Fs) + Fs, lambda = hazard, s = cutpoint)
  # }
  #
  # time_imp <- sapply(time, impute1)

  # Check: impute timed occur after landmark observed times
  if (any(time > time_imp)) {
    stop("Imputed times cannot precede the observed times!")
  }

  if (!is.null(maxtime)) {
    # If 'maxtime' is lower than observed time, censor the data
    min_time <- pmin(time_imp, maxtime)
    event <- as.numeric(time_imp == min_time)
    dat <- data.frame(time = min_time, event = event)
  } else {
    dat <- data.frame(time = time_imp, event = rep(1, length(time)))
  }

  return(dat)

}
