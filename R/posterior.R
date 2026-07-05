#' @title Posterior distribution of piecewise exponential constant hazard rates
#'
#' @description Using the Beta-Gamma conjugacy property, the posterior
#'   distribution of the piecewise hazard rates (\eqn{\lambda_j}, for `j = 1,
#'   ..., J`) is calculated and sampled from.
#'
#' @inheritParams survival_adapt
#' @inheritParams haz_to_prop
#' @param data data frame. Minimum requirements are 3 columns: event time
#'   (`time`), indicator of the event (`event`), and indicator for
#'   treatment assignment (`treatment`, coded `1` for treatment and
#'   `0` for control). Other columns can be included in the data frame and
#'   will be handled in the split.
#' @param empty_interval character. Policy for piecewise intervals with no
#'   exposed subjects in a treatment arm. `"propagate"` copies sufficient
#'   statistics from the nearest non-empty interval in the same arm;
#'   `"prior"` leaves the interval with zero exposure and zero events, making
#'   the posterior prior-driven; `"error"` stops with a clear message.
#'
#' @return An array of dimension 3. The first dimension is of length
#'   `N_mcmc`, the second dimension is of length \eqn{J} (one column for
#'   each hazard piece), and the third dimension is of length 2, with the first
#'   slice including posterior samples from `post_treatment`, and the
#'   second slice including posterior samples from `post_control`.
#'
#' @importFrom stats rgamma
#'
#' @noRd
posterior <- function(
  data,
  cutpoints,
  prior,
  N_mcmc,
  single_arm,
  empty_interval = "propagate"
) {
  empty_interval <- match.arg(empty_interval, c("propagate", "prior", "error"))
  n_intervals <- length(cutpoints)

  # Verify the expected treatment groups are actually present before
  # summarizing; when `treatment` is numeric, an absent group yields no summary
  # row at all, which would otherwise silently produce an all-NA posterior slice.
  if (sum(data$treatment == 1) == 0) {
    stop("No subjects in the treatment arm")
  }
  if (!single_arm && sum(data$treatment == 0) == 0) {
    stop("No subjects in the control arm")
  }

  data_summ <- posterior_sufficient_stats(data, cutpoints, single_arm)

  if (any(data_summ$n == 0)) {
    if (empty_interval == "error") {
      stop(
        "At least one treatment arm interval has zero subjects; set ",
        "'empty_interval' to 'propagate' or 'prior' to continue."
      )
    }

    if (empty_interval == "propagate") {
      data_summ <- propagate_empty_intervals(data_summ)
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

  if (!single_arm) {
    # If control patients present
    for (j in 1:n_intervals) {
      post[, j, 2] <- with(
        subset(data_summ, treatment == 0),
        rgamma(N_mcmc, prior[1] + tot_events[j], prior[2] + tot_time[j])
      )
    }
  }

  return(post)
}

propagate_empty_intervals <- function(data_summ) {
  for (treatment_value in unique(data_summ$treatment)) {
    treatment_rows <- which(data_summ$treatment == treatment_value)
    treatment_summ <- data_summ[treatment_rows, ]

    if (all(treatment_summ$n == 0)) {
      stop("No non-empty intervals with treatment = ", treatment_value)
    }

    # Walk the intervals in order; for an empty interval, carry forward the
    # last non-empty interval's data. If the first interval(s) are empty,
    # back-fill from the first non-empty interval instead.
    first_nonzero <- which(treatment_summ$n > 0)[1]
    for (k in seq_len(nrow(treatment_summ))) {
      if (treatment_summ$n[k] == 0) {
        source_k <- if (k > first_nonzero) k - 1 else first_nonzero
        warning(
          "Treatment value ",
          treatment_value,
          ", interval ",
          k,
          " has zero subjects; propagating data from interval ",
          source_k,
          " for posterior estimation.",
          call. = FALSE
        )
        treatment_summ[k, c("tot_time", "tot_events")] <-
          treatment_summ[source_k, c("tot_time", "tot_events")]
      }
    }

    data_summ[treatment_rows, c("tot_time", "tot_events")] <-
      treatment_summ[, c("tot_time", "tot_events")]
  }

  data_summ
}

posterior_sufficient_stats <- function(data, cutpoints, single_arm) {
  n_intervals <- length(cutpoints)
  interval_upper <- c(cutpoints[-1], Inf)
  treatment_values <- if (single_arm) 1 else c(0, 1)

  data_summ <- expand.grid(
    interval = seq_len(n_intervals),
    treatment = treatment_values
  )
  data_summ <- data_summ[c("treatment", "interval")]
  data_summ$n <- 0L
  data_summ$tot_time <- 0
  data_summ$tot_events <- 0

  for (treatment_value in treatment_values) {
    treatment_data <- data[data$treatment == treatment_value, , drop = FALSE]
    if (nrow(treatment_data) == 0) {
      next
    }

    for (j in seq_len(n_intervals)) {
      lower <- cutpoints[j]
      upper <- interval_upper[j]
      exposure <- pmax(0, pmin(treatment_data$time, upper) - lower)
      row <- data_summ$treatment == treatment_value & data_summ$interval == j

      data_summ$n[row] <- sum(exposure > 0)
      data_summ$tot_time[row] <- sum(exposure)
      data_summ$tot_events[row] <- sum(
        treatment_data$event == 1 &
          treatment_data$time > lower &
          treatment_data$time <= upper
      )
    }
  }

  data_summ$interval <- factor(
    data_summ$interval,
    levels = seq_len(n_intervals)
  )
  data_summ
}
