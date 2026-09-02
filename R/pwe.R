#' @title Simulate piecewise exponential time-to-event outcomes
#'
#' @description Simulates event times from a piecewise-exponential distribution,
#'   with optional administrative censoring at a fixed follow-up time.
#'
#' @param n A single non-negative integer giving the number of event times to
#'   simulate. The default is `1`; `n = 0` returns a zero-row data frame.
#' @param hazard A numeric vector of finite, non-negative event rates, with one
#'   value per interval defined by `cutpoints`. The default is `1`, giving a
#'   constant unit rate. If at least one outcome is requested and the final rate
#'   is zero, `maxtime` must be supplied so that subjects without an event can
#'   be administratively censored.
#' @param cutpoints `NULL` (the default), or a numeric vector of finite,
#'   positive, strictly increasing interior
#'   times at which the hazard rate changes. The number of hazard rates must be
#'   one greater than the number of cutpoints. Use `NULL` for a constant hazard.
#' @param maxtime `NULL` (the default), or a single finite, positive numeric
#'   administrative censoring time. When supplied, it must be later than every
#'   cutpoint.
#'
#' @details PWEALL represents the generating hazard with pieces closed on the
#'   left and open on the right. Because the event-time distribution is
#'   continuous, the value of
#'   the hazard at an isolated cutpoint does not alter the cumulative hazard,
#'   distribution, or generated samples. When realized event times are later
#'   assigned to analysis intervals, `goldilocks` follows the survival
#'   counting-process convention, open on the left and closed on the right, so
#'   an event exactly at a cutpoint belongs to the interval ending there. See
#'   [pwe_impute()] for the conditional sampling details.
#'
#' @return A data frame with one row per simulated subject and columns `time`,
#'   the event or censoring time, and `event`, coded `1` for an event and `0`
#'   for administrative censoring.
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
  validate_nonnegative_integer_scalar(n, "n")
  validate_cutpoints(cutpoints)
  validate_piecewise_hazard(hazard, cutpoints)
  validate_maxtime(maxtime, cutpoints)

  if (n == 0L) {
    return(data.frame(time = numeric(), event = numeric()))
  }

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
#' @description Draws an event time from a piecewise-exponential distribution
#'   conditional on a subject remaining event-free through the observed
#'   follow-up time, with optional administrative censoring.
#'
#' @inheritParams pwe_sim
#' @param time A required numeric vector of finite, non-negative event-free
#'   follow-up times for subjects who have not had an event. When `maxtime` is
#'   supplied, no value may exceed it. A zero-length vector returns
#'   a zero-row data frame. Values are not recycled against other arguments.
#' @param hazard A required numeric vector of finite, non-negative event rates,
#'   with one value per interval defined by `cutpoints`. If the final rate is
#'   zero, `maxtime` must be supplied.
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
#'   PWEALL represents the generating hazard with pieces closed on the left and
#'   open on the right. Its cumulative distribution is continuous at every
#'   cutpoint, so this endpoint choice does not affect imputation. For assigning
#'   realized event times to analysis intervals, `goldilocks` uses the survival
#'   counting-process convention, open on the left and closed on the right; an
#'   event exactly at a cutpoint belongs to the interval ending there.
#'
#' @return A data frame with one row per subject and columns `time`, the imputed
#'   event or censoring time, and `event`, coded `1` for an event and `0` for
#'   administrative censoring.
#' @export
#'
#' @examples
#' pwe_impute(time = c(3, 4, 5), hazard = c(0.002, 0.01), cutpoints = 12)
#' pwe_impute(time = c(3, 4, 5), hazard = c(0.002, 0.01), cutpoints = 12,
#'            maxtime = 36)
#' pwe_impute(time = 19.621870008, hazard = c(2.585924e-02, 3.685254e-09),
#'            cutpoints = 12, maxtime = 36)
pwe_impute <- function(time, hazard, cutpoints = NULL, maxtime = NULL) {
  validate_nonnegative_numeric_vector(time, "time", allow_empty = TRUE)
  validate_cutpoints(cutpoints)
  validate_piecewise_hazard(hazard, cutpoints)
  validate_maxtime(maxtime, cutpoints)

  if (length(time) == 0L) {
    return(data.frame(time = numeric(), event = numeric()))
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


#' @title Conditional event probability from a piecewise exponential model
#'
#' @description Calculates the probability of an event by a fixed endpoint,
#'   conditional on remaining event-free through an observed follow-up time.
#'
#' @inheritParams pwe_sim
#' @inheritParams survival_adapt
#' @param time A numeric vector of finite, non-negative event-free follow-up
#'   times.
#'
#' @return A numeric vector of conditional event probabilities.
#' @noRd
pwe_conditional_event_probability <- function(
  time,
  hazard,
  end_of_study,
  cutpoints = NULL
) {
  validate_cutpoints(cutpoints)
  validate_piecewise_hazard(hazard, cutpoints)
  validate_endpoint_time(end_of_study, cutpoints, "end_of_study")

  validate_nonnegative_numeric_vector(time, "time", allow_empty = TRUE)
  if (any(time > end_of_study)) {
    stop("'time' must not be later than 'end_of_study'")
  }
  if (length(time) == 0) {
    return(numeric())
  }

  interval_lower <- c(0, cutpoints)
  interval_upper <- c(cutpoints, Inf)
  remaining_exposure <- outer(
    time,
    seq_along(hazard),
    function(observed_time, interval) {
      pmax(
        0,
        pmin(end_of_study, interval_upper[interval]) -
          pmax(observed_time, interval_lower[interval])
      )
    }
  )
  remaining_cumulative_hazard <- drop(remaining_exposure %*% hazard)

  # P(T <= end_of_study | T > time) = 1 - exp(-(H(T*) - H(T))).
  # expm1 avoids cancellation when the remaining cumulative hazard is small.
  pmin(
    1,
    pmax(
      0,
      cumulative_hazard_to_probability(remaining_cumulative_hazard)
    )
  )
}


#' @title Calculate endpoint event probabilities from piecewise hazards
#'
#' @description Calculates the cumulative event probability at a fixed
#'   follow-up time for one or more sets of piecewise-constant hazard rates.
#'
#' @param hazard A required numeric matrix of finite, non-negative hazard rates.
#'   Rows represent parameter sets, such as posterior draws, and columns
#'   represent the intervals defined by `cutpoints`. The number of columns must
#'   equal `length(cutpoints) + 1`, and at least one row is required.
#' @param end_of_study A required single finite, positive numeric time at which
#'   the cumulative event
#'   probability is evaluated. It must be greater than every cutpoint.
#' @inheritParams pwe_sim
#' @inheritParams survival_adapt
#'
#' @details The cumulative probability depends on interval durations, so the
#'   value assigned to an isolated cutpoint has no effect. PWEALL represents
#'   its generating hazard with pieces closed on the left and open on the
#'   right. When `goldilocks` assigns realized event times to analysis
#'   intervals, it instead uses the survival counting-process convention,
#'   open on the left and closed on the right.
#'
#' @return A numeric vector of event probabilities in `[0, 1]`, with one value
#'   for each row of `hazard`.
#'
#' @export
ppwe <- function(hazard, end_of_study, cutpoints = NULL) {
  validate_cutpoints(cutpoints)
  validate_endpoint_time(end_of_study, cutpoints, "end_of_study")
  validate_hazard_matrix(hazard, cutpoints)

  duration <- endpoint_interval_widths(cutpoints, end_of_study)
  cumulative_hazard_to_probability(drop(hazard %*% duration))
}

#' Calculate analysis-interval widths through the endpoint horizon
#'
#' @description Returns the fixed model-interval durations used to integrate
#'   piecewise hazards through `end_of_study`. These widths describe the
#'   estimand horizon and do not depend on the maximum follow-up currently
#'   observed in the data.
#'
#' @param cutpoints `NULL`, or a numeric vector of finite, positive, strictly
#'   increasing analysis-model cutpoints.
#' @param end_of_study A single finite, positive numeric value giving the
#'   endpoint horizon.
#'
#' @return A numeric vector with one positive width per analysis interval.
#'
#' @keywords internal
#' @noRd
endpoint_interval_widths <- function(cutpoints, end_of_study) {
  diff(c(0, cutpoints, end_of_study))
}
