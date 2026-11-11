#' @title Analyze Bayesian survival sufficient statistics
#'
#' @description Applies the Bayesian piecewise-exponential final analysis
#'   directly to completed-data sufficient statistics. Predictive-imputation
#'   callers use this function after deriving the completed exposure and event
#'   totals, avoiding a second patient-level posterior setup.
#'
#' @param data_summ A data frame of completed-data exposure times and event
#'   counts by arm and piecewise interval, created by
#'   `posterior_sufficient_stats()`.
#' @inheritParams analyse_data
#'
#' @details This is the second posterior in the package's two-stage Bayesian
#'   procedure:
#'
#'   1. Draw hazards from the posterior conditional on the currently observed
#'      data.
#'   2. Use those hazards to impute a completed trial.
#'   3. Form a fresh posterior from the completed trial's sufficient statistics
#'      and the prespecified analysis prior (the caller supplies
#'      `prior_surv_final` as `prior_surv` for predictive and final analyses).
#'
#'   Step 3 deliberately does not update from the first posterior. The first
#'   posterior is used to generate the missing outcomes; treating it as the new
#'   prior would include the observed subjects both in that prior and again in
#'   `data_summ`.
#'
#' @return A list containing the posterior tail probability (`success`) and
#'   posterior mean event-probability effect (`effect`).
#'
#' @noRd
analyse_bayes_surv_sufficient_stats <- function(
  data_summ,
  cutpoints,
  end_of_study,
  prior_surv,
  N_mcmc,
  single_arm,
  alternative,
  h0,
  empty_interval = "prior"
) {
  validate_h0(h0, method = "bayes-surv", single_arm = single_arm)
  if (!alternative %in% c("greater", "less")) {
    stop(
      "Bayesian survival analysis requires 'alternative' to be ",
      "'greater' or 'less'"
    )
  }

  # Draw completed-data hazards directly from D_ak (events) and Y_ak
  # (exposure). posterior_from_sufficient_stats() combines these with the
  # original Gamma prior and preserves the configured empty-interval policy.
  post_lambda <- posterior_from_sufficient_stats(
    data_summ = data_summ,
    prior_surv = prior_surv,
    N_mcmc = N_mcmc,
    single_arm = single_arm,
    empty_interval = empty_interval
  )

  analyse_bayes_surv_posterior_draws(
    post_lambda = post_lambda,
    cutpoints = cutpoints,
    end_of_study = end_of_study,
    single_arm = single_arm,
    alternative = alternative,
    h0 = h0
  )
}

#' Analyze validated Bayesian survival sufficient statistics
#'
#' @description Applies the Bayesian survival analysis to sufficient statistics
#'   already validated and ordered by arm and piecewise interval.
#'
#' @inheritParams analyse_bayes_surv_sufficient_stats
#'
#' @return See `analyse_bayes_surv_sufficient_stats()`.
#'
#' @keywords internal
#' @noRd
analyse_prepared_bayes_surv <- function(
  data_summ,
  cutpoints,
  end_of_study,
  prior_surv,
  N_mcmc,
  single_arm,
  alternative,
  h0,
  empty_interval = "prior",
  interval_widths = NULL
) {
  if (
    !is.character(alternative) ||
      length(alternative) != 1L ||
      !alternative %in% c("greater", "less") ||
      !is.numeric(h0) ||
      length(h0) != 1L ||
      is.na(h0) ||
      !is.finite(h0)
  ) {
    stop(
      "Internal Bayesian-survival invariant failed: invalid hypothesis",
      call. = FALSE
    )
  }
  if (is.null(interval_widths)) {
    interval_widths <- endpoint_interval_widths(cutpoints, end_of_study)
  }
  valid_interval_widths <- is.numeric(interval_widths) &&
    length(interval_widths) == length(cutpoints) + 1L &&
    !anyNA(interval_widths) &&
    all(is.finite(interval_widths)) &&
    all(interval_widths > 0)
  if (!valid_interval_widths) {
    stop(
      "Internal Bayesian-survival invariant failed: invalid interval widths",
      call. = FALSE
    )
  }

  post_lambda <- posterior_from_prepared_stats(
    data_summ = data_summ,
    prior_surv = prior_surv,
    N_mcmc = N_mcmc,
    single_arm = single_arm,
    empty_interval = empty_interval
  )

  cumulative_hazards <- bayes_surv_cumulative_hazards(
    post_lambda = post_lambda,
    interval_widths = interval_widths,
    single_arm = single_arm
  )
  summarise_bayes_surv_hazards(
    cumulative_hazards = cumulative_hazards,
    alternative = alternative,
    h0 = h0
  )
}

#' Summarize Bayesian-survival posterior draws
#'
#' @param post_lambda A three-dimensional numeric array of posterior
#'   piecewise-hazard draws, with draws in rows, intervals in columns, and
#'   treatment followed by control in the third dimension.
#' @inheritParams analyse_bayes_surv_sufficient_stats
#'
#' @return See `analyse_bayes_surv_sufficient_stats()`.
#'
#' @keywords internal
#' @noRd
analyse_bayes_surv_posterior_draws <- function(
  post_lambda,
  cutpoints,
  end_of_study,
  single_arm,
  alternative,
  h0
) {
  # Use the same cumulative-hazard comparison as prepared predictive analyses.
  # Converting each arm to an event risk first can round both risks to one
  # and incorrectly turn a nonzero treatment difference into a tie.
  cumulative_hazards <- bayes_surv_cumulative_hazards(
    post_lambda = post_lambda,
    interval_widths = endpoint_interval_widths(cutpoints, end_of_study),
    single_arm = single_arm
  )

  summarise_bayes_surv_hazards(
    cumulative_hazards = cumulative_hazards,
    alternative = alternative,
    h0 = h0
  )
}

#' Integrate posterior hazard draws through the endpoint
#'
#' @description Retains cumulative hazards before converting them to event
#'   probabilities, so comparisons remain possible when those probabilities
#'   would round to one.
#'
#' @param interval_widths A numeric vector of positive analysis-interval
#'   durations through the endpoint horizon.
#' @inheritParams analyse_bayes_surv_posterior_draws
#'
#' @return A list with cumulative-hazard vectors named `treatment` and `control`.
#'   The control component is `NULL` for single-arm analyses.
#'
#' @keywords internal
#' @noRd
bayes_surv_cumulative_hazards <- function(
  post_lambda,
  interval_widths,
  single_arm
) {
  post_dimensions <- dim(post_lambda)
  valid_structure <- is.array(post_lambda) &&
    length(post_dimensions) == 3L &&
    post_dimensions[2L] == length(interval_widths) &&
    post_dimensions[3L] == 2L &&
    all(is.finite(interval_widths)) &&
    all(interval_widths > 0)
  if (!valid_structure) {
    stop(
      paste0(
        "Internal Bayesian-survival invariant failed: incompatible posterior ",
        "hazards and interval widths"
      ),
      call. = FALSE
    )
  }

  n_draws <- post_dimensions[1L]
  n_intervals <- post_dimensions[2L]
  treatment_hazard <- matrix(
    post_lambda[,, 1L],
    nrow = n_draws,
    ncol = n_intervals
  )
  treatment <- drop(treatment_hazard %*% interval_widths)
  control <- if (single_arm) {
    NULL
  } else {
    control_hazard <- matrix(
      post_lambda[,, 2L],
      nrow = n_draws,
      ncol = n_intervals
    )
    drop(control_hazard %*% interval_widths)
  }
  # Overflow must be reported as a failed calculation; two infinite hazards
  # do not provide enough numerical information to compare the treatment arms.
  if (
    any(!is.finite(treatment)) || any(treatment < 0) ||
      any(!is.finite(control)) || any(control < 0)
  ) {
    stop(
      "Bayesian survival analysis requires finite, non-negative cumulative hazards",
      call. = FALSE
    )
  }
  list(treatment = treatment, control = control)
}

#' Calculate a stable treatment-control event-risk difference
#'
#' @param treatment,control Numeric vectors of finite cumulative hazards.
#' @return The difference (`effect`), its `sign`, and its `log_absolute` value.
#' @noRd
bayes_surv_risk_difference <- function(treatment, control) {
  difference <- treatment - control
  gap <- abs(difference)
  # p_T - p_C = exp(-H_C) - exp(-H_T). Factoring out exp(-min(H_T,H_C))
  # avoids subtracting rounded event probabilities. Retain the log magnitude
  # even when the reported difference is too small to represent as a double.
  # For log(1 - exp(-gap)), expm1 preserves small gaps, while log1p preserves
  # the distance from one for large gaps. A genuine tie has log magnitude -Inf.
  log_absolute <- -pmin(treatment, control) + ifelse(
    gap <= log(2),
    log(-expm1(-gap)),
    log1p(-exp(-gap))
  )
  list(
    effect = sign(difference) * exp(log_absolute),
    sign = sign(difference),
    log_absolute = log_absolute
  )
}

#' Summarize Bayesian-survival posterior draws without rounded-risk comparisons
#'
#' @param cumulative_hazards A list returned by `bayes_surv_cumulative_hazards()`.
#' @inheritParams analyse_bayes_surv_posterior_draws
#'
#' @return See `analyse_bayes_surv_sufficient_stats()`.
#'
#' @keywords internal
#' @noRd
summarise_bayes_surv_hazards <- function(
  cumulative_hazards,
  alternative,
  h0
) {
  treatment <- cumulative_hazards$treatment
  control <- cumulative_hazards$control
  if (is.null(control)) {
    effect_draws <- cumulative_hazard_to_probability(treatment)
    # p_T < h0 is equivalent to H_T < -log(1-h0). Thresholds zero and one
    # map to zero and infinity, preserving strict comparisons at both ends.
    threshold <- probability_to_cumulative_hazard(h0)
    success_indicator <- if (alternative == "greater") {
      treatment > threshold
    } else {
      treatment < threshold
    }
  } else {
    difference <- bayes_surv_risk_difference(treatment, control)
    effect_draws <- difference$effect
    if (h0 == 0) {
      # At a zero risk-difference margin, cumulative-hazard ordering is exactly
      # equivalent, including under nonproportional piecewise hazards.
      success_indicator <- if (alternative == "greater") {
        treatment > control
      } else {
        treatment < control
      }
    } else if (abs(h0) == 1) {
      # With finite cumulative hazards the difference is strictly inside (-1,1),
      # even when its displayed value rounds to an endpoint of that interval.
      success_indicator <- rep.int(
        (h0 == 1) == (alternative == "less"), length(treatment)
      )
    } else if (h0 > 0) {
      # Nonzero margins remain risk differences. Compare signs first, then
      # log magnitudes; converting h0 to a hazard-difference margin is invalid.
      success_indicator <- if (alternative == "greater") {
        difference$sign > 0 & difference$log_absolute > log(h0)
      } else {
        difference$sign <= 0 | difference$log_absolute < log(h0)
      }
    } else {
      # Among negative differences, a larger magnitude means a smaller effect.
      success_indicator <- if (alternative == "greater") {
        difference$sign >= 0 | difference$log_absolute < log(-h0)
      } else {
        difference$sign < 0 & difference$log_absolute > log(-h0)
      }
    }
  }
  success <- mean(success_indicator)

  result <- list(
    success = success,
    effect = mean(effect_draws)
  )
  set_analysis_mc_counts(
    result,
    successes = sum(success_indicator),
    draws = length(success_indicator)
  )
}
