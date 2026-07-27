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
  prior_surv,
  N_mcmc,
  single_arm,
  empty_interval = "propagate"
) {
  empty_interval <- match.arg(empty_interval, c("propagate", "prior", "error"))
  n_intervals <- length(cutpoints) + 1L

  # Verify the expected treatment groups are actually present before
  # summarizing; when `treatment` is numeric, an absent group yields no summary
  # row at all, which would otherwise silently produce an all-NA posterior slice.
  if (sum(data$treatment == 1) == 0) {
    stop("No subjects in the treatment arm")
  }
  if (!single_arm && sum(data$treatment == 0) == 0) {
    stop("No subjects in the control arm")
  }

  # Keep the patient-level entry point as a thin compatibility wrapper. The
  # posterior itself depends only on the statistics below, so separating the
  # two operations lets repeated imputation paths calculate the statistics once
  # and pass them directly to posterior_from_sufficient_stats().
  data_summ <- posterior_sufficient_stats(data, cutpoints, single_arm)

  posterior_from_sufficient_stats(
    data_summ = data_summ,
    prior_surv = prior_surv,
    N_mcmc = N_mcmc,
    single_arm = single_arm,
    empty_interval = empty_interval
  )
}

#' @title Draw piecewise-exponential hazards from sufficient statistics
#'
#' @description Draws from the conjugate Gamma posterior using one row per
#'   treatment arm and piecewise interval. This is the statistics-based half of
#'   `posterior()`: callers that already have completed-data exposure times and
#'   event counts can avoid sending a patient-level data frame through the
#'   posterior calculation again.
#'
#' @param data_summ Data frame containing `treatment`, `interval`, `n`,
#'   `tot_time`, and `tot_events`. The `n` column is retained even though it is
#'   not part of the Gamma update because it distinguishes an empty interval
#'   from an interval having zero events. That distinction is required by
#'   `empty_interval`.
#' @inheritParams posterior
#'
#' @details For arm \eqn{a} and interval \eqn{k}, the posterior is
#'   \deqn{\lambda_{ak} \mid data \sim Gamma(\alpha + D_{ak},
#'   \beta + Y_{ak}),}
#'   where \eqn{D_{ak}} is `tot_events`, \eqn{Y_{ak}} is `tot_time`, and
#'   `prior_surv = c(alpha, beta)` uses the shape-rate parameterization.
#'
#'   Keeping this function separate is important for the package's two-stage
#'   Bayesian imputation procedure. A first posterior draw generates a completed
#'   data set; the sufficient statistics of that completed data set are then
#'   combined with the *original* prior here to form a fresh second posterior.
#'   Using the first posterior as the second-stage prior would count the observed
#'   data twice.
#'
#' @return See `posterior()`.
#'
#' @noRd
posterior_from_sufficient_stats <- function(
  data_summ,
  prior_surv,
  N_mcmc,
  single_arm,
  empty_interval = "propagate"
) {
  empty_interval <- match.arg(empty_interval, c("propagate", "prior", "error"))
  required_columns <- c(
    "treatment",
    "interval",
    "n",
    "tot_time",
    "tot_events"
  )
  missing_columns <- setdiff(required_columns, names(data_summ))
  if (length(missing_columns) > 0) {
    stop(
      "'data_summ' is missing required column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }

  # A complete sufficient-statistics table has exactly one row for every
  # expected arm/interval combination. Validate this here because malformed
  # tables could otherwise recycle values silently inside rgamma().
  treatment_values <- if (single_arm) 1 else c(0, 1)
  interval_values <- unique(as.character(data_summ$interval))
  n_intervals <- length(interval_values)
  prior_surv <- normalize_gamma_prior(
    prior_surv,
    n_intervals = n_intervals,
    name = "prior_surv"
  )
  expected_combinations <- expand.grid(
    treatment = treatment_values,
    interval = interval_values,
    stringsAsFactors = FALSE
  )
  actual_combinations <- paste(
    data_summ$treatment,
    as.character(data_summ$interval),
    sep = "\r"
  )
  expected_keys <- paste(
    expected_combinations$treatment,
    expected_combinations$interval,
    sep = "\r"
  )
  if (
    n_intervals == 0 ||
      nrow(data_summ) != nrow(expected_combinations) ||
      !setequal(actual_combinations, expected_keys) ||
      anyDuplicated(actual_combinations)
  ) {
    stop(
      "'data_summ' must contain exactly one row for every expected ",
      "treatment-arm and interval combination"
    )
  }

  statistic_columns <- c("n", "tot_time", "tot_events")
  statistic_values <- unlist(
    data_summ[statistic_columns],
    use.names = FALSE
  )
  if (
    anyNA(statistic_values) ||
      any(!is.finite(statistic_values)) ||
      any(statistic_values < 0)
  ) {
    stop("Posterior sufficient statistics must be finite and non-negative")
  }

  # `n` is a row count and `tot_events` is an event count. Requiring
  # integer-valued counts catches corrupt summaries before they reach the Gamma
  # sampler while still allowing the columns to use R's numeric storage mode.
  count_values <- unlist(
    data_summ[c("n", "tot_events")],
    use.names = FALSE
  )
  if (any(count_values != floor(count_values))) {
    stop("Posterior subject and event counts must be whole numbers")
  }

  # The statistics table deliberately contains zero-filled rows for expected
  # arms that are absent. Detect those here so the direct path fails with the
  # same diagnostic as posterior(data = ...), rather than eventually failing
  # inside empty-interval propagation.
  if (sum(data_summ$n[data_summ$treatment == 1]) == 0) {
    stop("No subjects in the treatment arm")
  }
  if (
    !single_arm &&
      sum(data_summ$n[data_summ$treatment == 0]) == 0
  ) {
    stop("No subjects in the control arm")
  }

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

  # The third dimension is always length two for compatibility with the
  # imputation code: slice 1 is treatment and slice 2 is control. In a
  # single-arm design the unused control slice deliberately remains NA.
  post <- array(NA_real_, dim = c(N_mcmc, n_intervals, 2))

  treatment_summ <- data_summ[data_summ$treatment == 1, , drop = FALSE]
  treatment_summ <- treatment_summ[
    match(interval_values, as.character(treatment_summ$interval)),
    ,
    drop = FALSE
  ]
  for (j in 1:n_intervals) {
    post[, j, 1] <- rgamma(
      N_mcmc,
      prior_surv[1, j] + treatment_summ$tot_events[j],
      prior_surv[2, j] + treatment_summ$tot_time[j]
    )
  }

  if (!single_arm) {
    control_summ <- data_summ[data_summ$treatment == 0, , drop = FALSE]
    control_summ <- control_summ[
      match(interval_values, as.character(control_summ$interval)),
      ,
      drop = FALSE
    ]
    for (j in 1:n_intervals) {
      post[, j, 2] <- rgamma(
        N_mcmc,
        prior_surv[1, j] + control_summ$tot_events[j],
        prior_surv[2, j] + control_summ$tot_time[j]
      )
    }
  }

  post
}

#' @title Propagate statistics for empty intervals
#'
#' @description Replaces empty piecewise intervals with nearby within-arm
#'   sufficient statistics under the historical propagation policy.
#'
#' @noRd
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

#' @title Calculate posterior sufficient statistics
#'
#' @description Aggregates exposure time and event counts by treatment arm and
#'   piecewise interval for conjugate Gamma posterior updates.
#'
#' @param rows Optional logical vector selecting the patient rows to summarize.
#'   This avoids constructing a temporary three-column analysis data frame when
#'   the caller already holds an imputed simulation data set.
#'
#' @noRd
posterior_sufficient_stats <- function(
  data,
  cutpoints,
  single_arm,
  rows = NULL
) {
  n_intervals <- length(cutpoints) + 1L
  interval_lower <- c(0, cutpoints)
  interval_upper <- c(cutpoints, Inf)
  treatment_values <- if (single_arm) 1 else c(0, 1)

  if (is.null(rows)) {
    rows <- rep.int(TRUE, nrow(data))
  } else if (
    !is.logical(rows) ||
      length(rows) != nrow(data) ||
      anyNA(rows)
  ) {
    stop("'rows' must be a non-missing logical vector with one value per row")
  }

  data_summ <- expand.grid(
    interval = seq_len(n_intervals),
    treatment = treatment_values
  )
  data_summ <- data_summ[c("treatment", "interval")]
  data_summ$n <- 0L
  data_summ$tot_time <- 0
  data_summ$tot_events <- 0

  for (treatment_value in treatment_values) {
    treatment_data <- data[
      rows & data$treatment == treatment_value,
      ,
      drop = FALSE
    ]
    if (nrow(treatment_data) == 0) {
      next
    }

    for (j in seq_len(n_intervals)) {
      lower <- interval_lower[j]
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
