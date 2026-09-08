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
#'      and the original `prior_surv`.
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

  effect_draws <- bayes_surv_effect_draws(
    post_lambda = post_lambda,
    interval_widths = interval_widths,
    single_arm = single_arm
  )
  summarise_bayes_surv_effect_draws(
    effect_draws = effect_draws,
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
  # haz_to_prop() evaluates all posterior draws together. For piecewise hazards
  # it delegates to the matrix-based ppwe() cumulative-hazard calculation,
  # avoiding a draw-by-draw probability loop.
  post_event_probability <- haz_to_prop(
    post = post_lambda,
    cutpoints = cutpoints,
    end_of_study = end_of_study,
    single_arm = single_arm
  )
  effect_draws <- post_event_probability$effect

  summarise_bayes_surv_effect_draws(
    effect_draws = effect_draws,
    alternative = alternative,
    h0 = h0
  )
}

#' Calculate Bayesian survival effects from posterior hazard draws
#'
#' @description Converts posterior piecewise-hazard draws to event probabilities
#'   at the endpoint horizon and, for a two-arm design, calculates the
#'   treatment-minus-control difference.
#'
#' @param interval_widths A numeric vector of positive analysis-interval
#'   durations through the endpoint horizon.
#' @inheritParams analyse_bayes_surv_posterior_draws
#'
#' @return A numeric vector containing one effect per posterior draw.
#'
#' @keywords internal
#' @noRd
bayes_surv_effect_draws <- function(
  post_lambda,
  interval_widths,
  single_arm
) {
  post_dimensions <- dim(post_lambda)
  valid_structure <- is.array(post_lambda) &&
    length(post_dimensions) == 3L &&
    post_dimensions[2L] == length(interval_widths) &&
    post_dimensions[3L] == 2L
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
  treatment_probability <- cumulative_hazard_to_probability(
    drop(treatment_hazard %*% interval_widths)
  )

  if (single_arm) {
    return(treatment_probability)
  }

  control_hazard <- matrix(
    post_lambda[,, 2L],
    nrow = n_draws,
    ncol = n_intervals
  )
  control_probability <- cumulative_hazard_to_probability(
    drop(control_hazard %*% interval_widths)
  )
  treatment_probability - control_probability
}

#' Summarize Bayesian-survival effect draws
#'
#' @param effect_draws A numeric vector of posterior event-probability effects.
#' @inheritParams analyse_bayes_surv_posterior_draws
#'
#' @return See `analyse_bayes_surv_sufficient_stats()`.
#'
#' @keywords internal
#' @noRd
summarise_bayes_surv_effect_draws <- function(
  effect_draws,
  alternative,
  h0
) {
  success_indicator <- if (alternative == "greater") {
    effect_draws > h0
  } else {
    effect_draws < h0
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
