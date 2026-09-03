#' @title Sample posterior piecewise-exponential hazard rates
#'
#' @description Updates independent Gamma priors with the observed event counts
#'   and exposure times, then samples the posterior distribution of each
#'   piecewise-constant hazard rate.
#'
#' @inheritParams survival_adapt
#' @inheritParams haz_to_prop
#' @param data A data frame with one row per subject and at least three columns:
#'   follow-up time (`time`), event indicator (`event`), and
#'   treatment assignment (`treatment`, coded `1` for treatment and
#'   `0` for control). Additional columns are ignored.
#' @param empty_interval A single character string specifying how to handle a
#'   piecewise interval with no
#'   exposed subjects in a treatment arm. `"propagate"` copies sufficient
#'   statistics from the nearest non-empty interval in the same arm;
#'   `"prior"` (the default) leaves the interval with zero exposure and zero
#'   events, making
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
  empty_interval = "prior"
) {
  empty_interval <- match.arg(empty_interval, c("prior", "propagate", "error"))
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
#' @description Draws from the conjugate Gamma posterior using one row of
#'   sufficient statistics per treatment arm and piecewise interval.
#'
#' @param data_summ A data frame containing `treatment`, `interval`, `n`,
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
  empty_interval = "prior"
) {
  empty_interval <- match.arg(empty_interval, c("prior", "propagate", "error"))
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
  interval_values <- if (is.factor(data_summ$interval)) {
    levels(data_summ$interval)
  } else {
    sort(unique(data_summ$interval))
  }
  interval_values <- as.character(interval_values)
  n_intervals <- length(interval_values)
  prior_surv <- normalize_gamma_prior(
    prior_surv,
    n_intervals = n_intervals,
    single_arm = single_arm,
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

  data_summ <- resolve_empty_posterior_intervals(
    data_summ,
    empty_interval = empty_interval
  )
  canonical_combinations <- expand.grid(
    interval = interval_values,
    treatment = treatment_values,
    stringsAsFactors = FALSE
  )
  canonical_keys <- paste(
    canonical_combinations$treatment,
    canonical_combinations$interval,
    sep = "\r"
  )
  data_summ <- data_summ[
    match(canonical_keys, actual_combinations),
    ,
    drop = FALSE
  ]

  draw_gamma_posterior_kernel(
    data_summ = data_summ,
    prior_surv = prior_surv,
    N_mcmc = N_mcmc,
    single_arm = single_arm
  )
}

#' Draw posterior hazards from validated sufficient statistics
#'
#' @description Draws posterior hazards from sufficient statistics already
#'   validated and ordered by arm and interval, then applies the selected
#'   empty-interval rule before sampling.
#'
#' @inheritParams posterior_from_sufficient_stats
#'
#' @return See `posterior()`.
#'
#' @keywords internal
#' @noRd
posterior_from_sufficient_stats_kernel <- function(
  data_summ,
  prior_surv,
  N_mcmc,
  single_arm,
  empty_interval = "prior"
) {
  if (
    !is.logical(single_arm) ||
      length(single_arm) != 1L ||
      is.na(single_arm)
  ) {
    stop(
      "Internal posterior invariant failed: invalid single-arm indicator",
      call. = FALSE
    )
  }
  if (
    !is.numeric(N_mcmc) ||
      length(N_mcmc) != 1L ||
      is.na(N_mcmc) ||
      !is.finite(N_mcmc) ||
      N_mcmc < 1L ||
      N_mcmc != floor(N_mcmc)
  ) {
    stop(
      "Internal posterior invariant failed: invalid draw count",
      call. = FALSE
    )
  }
  if (
    !is.character(empty_interval) ||
      length(empty_interval) != 1L ||
      !empty_interval %in% c("prior", "propagate", "error")
  ) {
    stop(
      "Internal posterior invariant failed: invalid empty-interval policy",
      call. = FALSE
    )
  }

  expected_arms <- if (single_arm) "treatment" else c("control", "treatment")
  valid_prior <- is.array(prior_surv) &&
    length(dim(prior_surv)) == 3L &&
    dim(prior_surv)[1L] == 2L &&
    dim(prior_surv)[3L] == length(expected_arms) &&
    identical(dimnames(prior_surv)$parameter, c("shape", "rate")) &&
    identical(dimnames(prior_surv)$arm, expected_arms) &&
    is.numeric(prior_surv) &&
    !anyNA(prior_surv) &&
    all(is.finite(prior_surv)) &&
    all(prior_surv > 0)
  if (!valid_prior) {
    stop(
      "Internal posterior invariant failed: non-canonical survival prior",
      call. = FALSE
    )
  }

  n_intervals <- dim(prior_surv)[2L]
  required_columns <- c(
    "treatment",
    "interval",
    "n",
    "tot_time",
    "tot_events"
  )
  expected_treatment <- if (single_arm) {
    rep(1, n_intervals)
  } else {
    rep(c(0, 1), each = n_intervals)
  }
  expected_interval <- rep(seq_len(n_intervals), length(expected_arms))
  valid_structure <- is.data.frame(data_summ) &&
    all(required_columns %in% names(data_summ)) &&
    nrow(data_summ) == length(expected_treatment) &&
    !anyNA(data_summ$treatment) &&
    all(data_summ$treatment == expected_treatment) &&
    !anyNA(data_summ$interval) &&
    all(as.character(data_summ$interval) == as.character(expected_interval))
  if (!valid_structure) {
    stop(
      paste0(
        "Internal posterior invariant failed: non-canonical sufficient ",
        "statistics"
      ),
      call. = FALSE
    )
  }

  statistic_values <- unlist(
    data_summ[c("n", "tot_time", "tot_events")],
    use.names = FALSE
  )
  count_values <- unlist(
    data_summ[c("n", "tot_events")],
    use.names = FALSE
  )
  if (
    anyNA(statistic_values) ||
      any(!is.finite(statistic_values)) ||
      any(statistic_values < 0) ||
      any(count_values != floor(count_values))
  ) {
    stop(
      "Internal posterior invariant failed: invalid sufficient statistics",
      call. = FALSE
    )
  }
  if (sum(data_summ$n[data_summ$treatment == 1]) == 0) {
    stop(
      "Internal posterior invariant failed: no subjects in treatment arm",
      call. = FALSE
    )
  }
  if (!single_arm && sum(data_summ$n[data_summ$treatment == 0]) == 0) {
    stop(
      "Internal posterior invariant failed: no subjects in control arm",
      call. = FALSE
    )
  }

  data_summ <- resolve_empty_posterior_intervals(
    data_summ,
    empty_interval = empty_interval
  )
  draw_gamma_posterior_kernel(
    data_summ = data_summ,
    prior_surv = prior_surv,
    N_mcmc = N_mcmc,
    single_arm = single_arm
  )
}

#' Apply the configured empty-interval policy
#'
#' @inheritParams posterior_from_sufficient_stats
#'
#' @return The sufficient-statistics data frame after any propagation.
#'
#' @keywords internal
#' @noRd
resolve_empty_posterior_intervals <- function(data_summ, empty_interval) {
  if (!any(data_summ$n == 0)) {
    return(data_summ)
  }

  empty_details <- data_summ[
    data_summ$n == 0,
    c("treatment", "interval"),
    drop = FALSE
  ]
  empty_details$interval <- as.character(empty_details$interval)
  signalCondition(structure(
    list(
      message = paste0(
        "Empty piecewise-exponential interval(s) handled with policy '",
        empty_interval,
        "'"
      ),
      call = NULL,
      policy = empty_interval,
      details = empty_details
    ),
    class = c("goldilocks_empty_interval", "condition")
  ))

  if (empty_interval == "error") {
    stop(
      "At least one treatment arm interval has zero subjects; set ",
      "'empty_interval' to 'propagate' or 'prior' to continue."
    )
  }
  if (empty_interval == "propagate") {
    return(propagate_empty_intervals(data_summ))
  }
  data_summ
}

#' Sample Gamma posterior hazards from canonical sufficient statistics
#'
#' @inheritParams posterior_from_sufficient_stats
#'
#' @return See `posterior()`.
#'
#' @keywords internal
#' @noRd
draw_gamma_posterior_kernel <- function(
  data_summ,
  prior_surv,
  N_mcmc,
  single_arm
) {
  n_intervals <- dim(prior_surv)[2L]
  treatment_rows <- if (single_arm) {
    seq_len(n_intervals)
  } else {
    n_intervals + seq_len(n_intervals)
  }
  control_rows <- seq_len(n_intervals)

  # The third dimension remains length two for compatibility with imputation:
  # slice 1 is treatment and slice 2 is control. The unused single-arm control
  # slice deliberately remains NA.
  post <- array(NA_real_, dim = c(N_mcmc, n_intervals, 2L))
  for (j in seq_len(n_intervals)) {
    row <- treatment_rows[j]
    post[, j, 1L] <- rgamma(
      N_mcmc,
      prior_surv["shape", j, "treatment"] + data_summ$tot_events[row],
      prior_surv["rate", j, "treatment"] + data_summ$tot_time[row]
    )
  }

  if (!single_arm) {
    for (j in seq_len(n_intervals)) {
      row <- control_rows[j]
      post[, j, 2L] <- rgamma(
        N_mcmc,
        prior_surv["shape", j, "control"] + data_summ$tot_events[row],
        prior_surv["rate", j, "control"] + data_summ$tot_time[row]
      )
    }
  }
  post
}

#' Summarize resolved Gamma priors by arm and interval
#'
#' @param prior_surv A numeric vector, matrix, or named list specifying a Gamma
#'   prior in one of the forms described for [survival_adapt()].
#' @param cutpoints `NULL`, or a numeric vector of piecewise-exponential
#'   cutpoints.
#' @param end_of_study A single finite, positive numeric value giving the
#'   subject-level follow-up horizon.
#' @param single_arm A single logical value indicating whether the design is
#'   single-arm.
#' @param stage `NULL`, or a single character string labelling the interim or
#'   final analysis stage. The default is `NULL`.
#'
#' @return A tidy data frame containing resolved Gamma parameters and moments.
#'
#' @keywords internal
#' @noRd
gamma_prior_diagnostics <- function(
  prior_surv,
  cutpoints,
  end_of_study,
  single_arm,
  stage = NULL
) {
  n_intervals <- length(cutpoints) + 1L
  prior_surv <- normalize_gamma_prior(
    prior_surv,
    n_intervals = n_intervals,
    single_arm = single_arm,
    name = "prior_surv"
  )
  arm_names <- dimnames(prior_surv)$arm
  interval_start <- c(0, cutpoints)
  interval_end <- c(cutpoints, end_of_study)
  rows <- lapply(arm_names, function(arm) {
    shape <- prior_surv["shape", , arm]
    rate <- prior_surv["rate", , arm]
    data.frame(
      arm = arm,
      interval = seq_len(n_intervals),
      interval_start = interval_start,
      interval_end = interval_end,
      shape = unname(shape),
      rate = unname(rate),
      mean_hazard = unname(shape / rate),
      sd_hazard = unname(sqrt(shape) / rate),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  if (!is.null(stage)) {
    out$stage <- stage
    out <- out[c("stage", setdiff(names(out), "stage"))]
  }
  out
}

#' Summarize conjugate Gamma posterior parameters
#'
#' @description Combines a resolved arm-specific prior with observed sufficient
#'   statistics. Under the legacy propagation policy, both the observed and
#'   effective statistics are retained so the resulting posterior remains
#'   auditable.
#'
#' @inheritParams gamma_prior_diagnostics
#' @param data_summ A data frame of observed piecewise-exponential sufficient
#'   statistics, with one row per arm and interval.
#' @param empty_interval A single character string specifying the empty-interval
#'   rule used by the posterior update.
#'
#' @return A tidy data frame containing prior parameters, observed and effective
#'   sufficient statistics, and posterior parameters and moments.
#'
#' @keywords internal
#' @noRd
gamma_posterior_diagnostics <- function(
  data_summ,
  prior_surv,
  cutpoints,
  end_of_study,
  single_arm,
  empty_interval
) {
  empty_interval <- match.arg(
    empty_interval,
    c("prior", "propagate", "error")
  )
  prior <- gamma_prior_diagnostics(
    prior_surv = prior_surv,
    cutpoints = cutpoints,
    end_of_study = end_of_study,
    single_arm = single_arm
  )
  effective_summ <- data_summ
  if (empty_interval == "propagate" && any(effective_summ$n == 0L)) {
    # The posterior update has already emitted the user-facing fallback
    # warning. Reproduce its effective statistics without warning twice.
    effective_summ <- suppressWarnings(
      propagate_empty_intervals(effective_summ)
    )
  }

  treatment_value <- ifelse(prior$arm == "treatment", 1, 0)
  observed_key <- paste(
    data_summ$treatment,
    as.character(data_summ$interval),
    sep = "\r"
  )
  diagnostic_key <- paste(treatment_value, prior$interval, sep = "\r")
  observed_rows <- match(diagnostic_key, observed_key)
  effective_key <- paste(
    effective_summ$treatment,
    as.character(effective_summ$interval),
    sep = "\r"
  )
  effective_rows <- match(diagnostic_key, effective_key)

  observed_n <- data_summ$n[observed_rows]
  observed_exposure <- data_summ$tot_time[observed_rows]
  observed_events <- data_summ$tot_events[observed_rows]
  effective_exposure <- effective_summ$tot_time[effective_rows]
  effective_events <- effective_summ$tot_events[effective_rows]
  posterior_shape <- prior$shape + effective_events
  posterior_rate <- prior$rate + effective_exposure

  data.frame(
    arm = prior$arm,
    interval = prior$interval,
    interval_start = prior$interval_start,
    interval_end = prior$interval_end,
    exposed_subjects = observed_n,
    observed_exposure = observed_exposure,
    observed_events = observed_events,
    empty_interval = observed_n == 0L,
    empty_interval_policy = ifelse(
      observed_n == 0L,
      empty_interval,
      "observed"
    ),
    effective_exposure = effective_exposure,
    effective_events = effective_events,
    prior_shape = prior$shape,
    prior_rate = prior$rate,
    posterior_shape = posterior_shape,
    posterior_rate = posterior_rate,
    posterior_mean_hazard = posterior_shape / posterior_rate,
    posterior_sd_hazard = sqrt(posterior_shape) / posterior_rate,
    stringsAsFactors = FALSE
  )
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
#' @param rows `NULL` (the default), or a logical vector with one value per row
#'   of `data`, selecting the subjects to include in the sufficient statistics.
#'
#' @noRd
posterior_sufficient_stats <- function(
  data,
  cutpoints,
  single_arm,
  rows = NULL
) {
  required_columns <- c("time", "event", "treatment")
  missing_columns <- setdiff(required_columns, names(data))
  if (!is.data.frame(data) || length(missing_columns) > 0L) {
    stop(
      "'data' must be a data frame containing columns: time, event, treatment"
    )
  }
  validate_cutpoints(cutpoints)
  validate_nonnegative_numeric_vector(data$time, "data$time")
  validate_binary_vector(data$event, "data$event")
  validate_binary_vector(data$treatment, "data$treatment")

  if (is.null(rows)) {
    rows <- rep.int(TRUE, nrow(data))
  } else if (
    !is.logical(rows) ||
      length(rows) != nrow(data) ||
      anyNA(rows)
  ) {
    stop("'rows' must be a non-missing logical vector with one value per row")
  }

  posterior_sufficient_stats_kernel(
    time = data$time,
    event = data$event,
    treatment = data$treatment,
    cutpoints = cutpoints,
    single_arm = single_arm,
    rows = rows
  )
}

#' Calculate posterior sufficient statistics from outcome vectors
#'
#' @description Aggregates exposure time and event counts by treatment arm and
#'   piecewise interval after the input vectors have already been validated.
#'
#' @param time A numeric vector of follow-up times.
#' @param event A binary vector of event indicators.
#' @param treatment A binary vector of treatment assignments.
#' @param cutpoints A numeric vector of interior piecewise-interval boundaries,
#'   or `NULL` for a constant hazard.
#' @param single_arm A single logical value indicating a one-arm design.
#' @param rows A logical vector selecting the subjects to include.
#'
#' @return A data frame containing subject counts, total exposure, and event
#'   counts by treatment arm and piecewise interval.
#'
#' @keywords internal
#' @noRd
posterior_sufficient_stats_kernel <- function(
  time,
  event,
  treatment,
  cutpoints,
  single_arm,
  rows = rep.int(TRUE, length(time))
) {
  n_intervals <- length(cutpoints) + 1L
  interval_lower <- c(0, cutpoints)
  interval_upper <- c(cutpoints, Inf)
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
    treatment_rows <- rows & treatment == treatment_value
    if (!any(treatment_rows)) {
      next
    }
    treatment_time <- time[treatment_rows]
    treatment_event <- event[treatment_rows]
    # Match Surv(start, stop, event) and survSplit(): a realized event at a
    # cutpoint belongs to the interval ending at that cutpoint.
    event_interval <- right_closed_interval_index(
      treatment_time,
      interval_lower
    )

    for (j in seq_len(n_intervals)) {
      lower <- interval_lower[j]
      upper <- interval_upper[j]
      exposure <- pmax(0, pmin(treatment_time, upper) - lower)
      row <- data_summ$treatment == treatment_value & data_summ$interval == j

      data_summ$n[row] <- sum(exposure > 0)
      data_summ$tot_time[row] <- sum(exposure)
      data_summ$tot_events[row] <- sum(
        treatment_event == 1 &
          event_interval == j
      )
    }
  }

  data_summ$interval <- factor(
    data_summ$interval,
    levels = seq_len(n_intervals)
  )
  data_summ
}
