#' @title Posterior distribution of piecewise exponential constant hazard rates
#'
#' @description Using the Beta-Gamma conjugacy property, the posterior
#'   distribution of the piecewise hazard rates (\eqn{\lambda_j}, for \code{j=1,
#'   \dots, J}) is calculated and sampled from.
#'
#' @inheritParams survival_adapt
#' @inheritParams haz_to_prop
#' @param data data frame. Minimum requirements are 3 columns: event time
#'   (\code{time}), indicator of the event (\code{event}), and indicator for
#'   treatment arm (\code{treatment}). Other columns can be included in the data
#'   frame and will be handled in the split.
#'
#' @return An array of dimension 3. The first dimension is of length
#'   \code{N_mcmc}, the second dimension is of length \eqn{J} (one column for
#'   each hazard piece), and the third dimension is of length 2, with the first
#'   slice including posterior samples from \code{post_treatment}, and the
#'   second slice including posterior samples from \code{post_control}.
#'
#' @importFrom stats rgamma
#'
#' @noRd
posterior <- function(data, cutpoints, prior, N_mcmc, single_arm) {

  n_intervals <- length(cutpoints)

  # Verify the expected treatment arms are actually present before summarising;
  # when `treatment` is numeric, an absent arm yields no summary row at all,
  # which would otherwise silently produce an all-NA posterior for that arm via
  # rgamma(N_mcmc, ..., NA).
  if (sum(data$treatment == 1) == 0) {
    stop("No subjects in the treatment arm")
  }
  if (!single_arm && sum(data$treatment == 0) == 0) {
    stop("No subjects in the control arm")
  }

  data_summ <- posterior_sufficient_stats(data, cutpoints, single_arm)

  # If a time-interval has zero subjects at a given interim analysis, it will
  # mean that there is zero exposure time and events for that stratum. To avoid
  # unrealistic hazard parameter estimates, we propagate the exposure time and
  # event counts from the nearest non-zero stratum *within the same treatment
  # arm*. This is equivalent to independent draws from the Gamma posterior of
  # that stratum. Propagation is done per arm because data_summ stacks the arms
  # (all control intervals, then all treatment intervals); propagating across
  # the flat row order would otherwise leak one arm's data into the other.
  if (any(data_summ$n == 0)) {
    for (arm in unique(data_summ$treatment)) {
      arm_rows <- which(data_summ$treatment == arm)
      arm_summ <- data_summ[arm_rows, ]

      if (all(arm_summ$n == 0)) {
        # Defensive: with the top-level arm-presence check above this should be
        # unreachable, but guard against an arm that appears only as empty
        # factor levels.
        stop("No subjects in arm ", arm)
      }

      # Walk the intervals in order; for an empty interval, carry forward the
      # last non-empty interval's data. If the first interval(s) are empty,
      # back-fill from the first non-empty interval instead.
      first_nonzero <- which(arm_summ$n > 0)[1]
      for (k in seq_len(nrow(arm_summ))) {
        if (arm_summ$n[k] == 0) {
          source_k <- if (k > first_nonzero) k - 1 else first_nonzero
          warning(
            "Treatment arm ", arm, ", interval ", k,
            " has zero subjects; propagating data from interval ", source_k,
            " for posterior estimation.",
            call. = FALSE
          )
          arm_summ[k, c("tot_time", "tot_events")] <-
            arm_summ[source_k, c("tot_time", "tot_events")]
        }
      }

      data_summ[arm_rows, c("tot_time", "tot_events")] <-
        arm_summ[, c("tot_time", "tot_events")]
    }
  }

  post <- array(dim = c(N_mcmc, n_intervals, 2))
  post_treatment <- matrix(nrow = N_mcmc, ncol = n_intervals)
  post_control <- matrix(nrow = N_mcmc, ncol = n_intervals)

  for (j in 1:n_intervals) {
    post[, j, 1] <- with(
      subset(data_summ, treatment == 1),
      rgamma(N_mcmc, prior[1] + tot_events[j], prior[2] + tot_time[j])
    )
  }

  if (!single_arm) { # If control patients present
    for (j in 1:n_intervals) {
      post[, j, 2] <- with(
        subset(data_summ, treatment == 0),
        rgamma(N_mcmc, prior[1] + tot_events[j], prior[2] + tot_time[j]))
    }
  }

  return(post)

}

posterior_sufficient_stats <- function(data, cutpoints, single_arm) {
  n_intervals <- length(cutpoints)
  interval_upper <- c(cutpoints[-1], Inf)
  arms <- if (single_arm) 1 else c(0, 1)

  data_summ <- expand.grid(
    interval = seq_len(n_intervals),
    treatment = arms
  )
  data_summ <- data_summ[c("treatment", "interval")]
  data_summ$n <- 0L
  data_summ$tot_time <- 0
  data_summ$tot_events <- 0

  for (arm in arms) {
    arm_data <- data[data$treatment == arm, , drop = FALSE]
    if (nrow(arm_data) == 0) {
      next
    }

    for (j in seq_len(n_intervals)) {
      lower <- cutpoints[j]
      upper <- interval_upper[j]
      exposure <- pmax(0, pmin(arm_data$time, upper) - lower)
      row <- data_summ$treatment == arm & data_summ$interval == j

      data_summ$n[row] <- sum(exposure > 0)
      data_summ$tot_time[row] <- sum(exposure)
      data_summ$tot_events[row] <- sum(
        arm_data$event == 1 &
          arm_data$time > lower &
          arm_data$time <= upper
      )
    }
  }

  data_summ$interval <- factor(data_summ$interval, levels = seq_len(n_intervals))
  data_summ
}
