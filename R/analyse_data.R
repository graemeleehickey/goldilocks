#' @title Apply the prespecified analysis to completed trial data
#'
#' @description Applies the selected Bayesian or frequentist analysis to an
#'   observed or imputed completed data set and returns its success measure and
#'   treatment-effect estimate.
#'
#' @inheritParams survival_adapt
#' @inheritParams sim_comp_data
#' @inheritParams haz_to_prop
#' @param data A data frame with one row per subject and columns `time`, `event`,
#'   and `treatment`. It contains the observed or imputed outcomes to be analyzed
#'   using the prespecified method.
#'
#' @return A list with two elements:
#'
#'   - `success`: Analysis-specific success score:
#'     - if `method = "bayes-surv"`, the posterior probability that the treatment
#'       effect is greater than `h0` when `alternative = "greater"`, or less
#'       than `h0` when `alternative = "less"`;
#'     - if `method = "logrank"`, 1 minus the log-rank test *P*-value, using a
#'       two-sided *P*-value when `alternative = "two.sided"` and a one-sided
#'       *P*-value otherwise;
#'     - if `method = "cox"`, 1 minus the Cox Wald test *P*-value for the
#'       estimated log hazard ratio compared with `h0`, using a two-sided
#'       *P*-value when `alternative = "two.sided"` and a one-sided *P*-value
#'       otherwise;
#'     - if `method = "bayes-bin"`, the posterior probability that the binary
#'       event proportion (single-arm) or treatment-control difference in
#'       binary event proportions (two-arm) is greater than `h0` when
#'       `alternative = "greater"`, or less than `h0` when
#'       `alternative = "less"`;
#'     - if `method = "riskdiff"`, 1 minus the Wald-test *P*-value for the
#'       treatment-control difference in binary event proportions compared
#'       with `h0`.
#'   - `effect`: Posterior mean effect for `method = "bayes-surv"` or
#'     `method = "bayes-bin"`, the estimated log hazard ratio for `method =
#'     "cox"`, the estimated treatment-control event-proportion difference for
#'     `method = "riskdiff"`, or `NA` for `method = "logrank"`.
#'
#' @importFrom stats dbeta integrate pbeta pnorm rbeta
#' @import Rcpp
#' @import survival
#' @useDynLib goldilocks, .registration = TRUE
#'
#' @noRd
analyse_data <- function(
  data,
  cutpoints,
  end_of_study,
  prior_surv,
  N_mcmc,
  single_arm,
  method,
  alternative,
  h0,
  prior_bin = c(1, 1),
  bin_method = "mc",
  empty_interval = "prior"
) {
  validate_h0(h0, method, single_arm)

  ####################################################
  ### Bayesian survival estimate test
  ####################################################

  # CIF_trt(T) - CIF_con(T) for two-armed trial

  if (method == "bayes-surv") {
    # The patient-level entry point computes the sufficient statistics once,
    # then delegates to the same fast path used after predictive imputation.
    # Keeping a single implementation prevents the ordinary and optimized
    # analysis paths from drifting apart statistically.
    data_summ <- posterior_sufficient_stats(
      data = data,
      cutpoints = cutpoints,
      single_arm = single_arm
    )
    return(analyse_bayes_surv_sufficient_stats(
      data_summ = data_summ,
      cutpoints = cutpoints,
      end_of_study = end_of_study,
      prior_surv = prior_surv,
      N_mcmc = N_mcmc,
      single_arm = single_arm,
      alternative = alternative,
      h0 = h0,
      empty_interval = empty_interval
    ))
  }

  ####################################################
  ### Bayesian binomial test
  ####################################################

  if (method == "bayes-bin") {
    assert_complete_binary_outcomes(data, end_of_study, "Bayesian binomial")
    bin_res <- bayes_binomial_test(
      data = data,
      single_arm = single_arm,
      alternative = alternative,
      h0 = h0,
      prior_bin = prior_bin,
      bin_method = bin_method,
      N_mcmc = N_mcmc
    )
    success <- bin_res$success
    effect <- bin_res$effect
  }

  ####################################################
  ### Log-rank test
  ####################################################

  if (method == "logrank") {
    analysis <- analyse_logrank_vectors(
      time = data$time,
      event = data$event,
      treatment = data$treatment,
      alternative = alternative
    )
    success <- analysis$success
    effect <- analysis$effect
  }

  ####################################################
  ### Cox regression test
  ####################################################

  if (method == "cox") {
    analysis <- analyse_cox_vectors(
      time = data$time,
      event = data$event,
      treatment = data$treatment,
      alternative = alternative,
      h0 = h0
    )
    success <- analysis$success
    effect <- analysis$effect
  }

  ####################################################
  ### Frequentist risk-difference test
  ####################################################

  if (method == "riskdiff") {
    fit_riskdiff <- risk_difference_wald_test_checked(
      data = data,
      end_of_study = end_of_study,
      alternative = alternative,
      h0 = h0
    )
    success <- fit_riskdiff$success
    effect <- fit_riskdiff$estimate
  }

  result <- list(
    "success" = success,
    "effect" = effect
  )
  if (method == "bayes-bin") {
    attr(result, "mc_counts") <- attr(bin_res, "mc_counts", exact = TRUE)
  }
  return(result)
}

#' Analyze completed survival outcomes from vectors
#'
#' @description Applies the selected completed-data survival analysis directly
#'   to follow-up, event, and treatment vectors. Predictive calculations use
#'   this function after inserting one imputation into the observed outcomes.
#'
#' @param outcome A list containing `time`, `event`, and `treatment` vectors.
#' @param interval_widths A numeric vector giving the piecewise-interval widths
#'   through `end_of_study` for a Bayesian survival analysis.
#' @inheritParams analyse_data
#'
#' @return A list containing the analysis-specific success score and treatment-
#'   effect estimate.
#'
#' @keywords internal
#' @noRd
analyse_completed_survival_vectors <- function(
  outcome,
  cutpoints,
  end_of_study,
  interval_widths,
  prior_surv,
  N_mcmc,
  single_arm,
  method,
  alternative,
  h0,
  empty_interval
) {
  if (method == "logrank") {
    return(analyse_logrank_vectors(
      time = outcome$time,
      event = outcome$event,
      treatment = outcome$treatment,
      alternative = alternative
    ))
  }
  if (method == "cox") {
    return(analyse_cox_vectors(
      time = outcome$time,
      event = outcome$event,
      treatment = outcome$treatment,
      alternative = alternative,
      h0 = h0
    ))
  }
  if (method == "bayes-surv") {
    data_summ <- posterior_sufficient_stats_kernel(
      time = outcome$time,
      event = outcome$event,
      treatment = outcome$treatment,
      cutpoints = cutpoints,
      single_arm = single_arm
    )
    return(analyse_bayes_surv_sufficient_stats_kernel(
      data_summ = data_summ,
      cutpoints = cutpoints,
      end_of_study = end_of_study,
      prior_surv = prior_surv,
      N_mcmc = N_mcmc,
      single_arm = single_arm,
      alternative = alternative,
      h0 = h0,
      empty_interval = empty_interval,
      interval_widths = interval_widths
    ))
  }

  stop(
    "Internal predictive-analysis invariant failed: unsupported survival method",
    call. = FALSE
  )
}

#' Calculate a log-rank result from outcome vectors
#'
#' @param time A numeric vector of follow-up times.
#' @param event A binary vector of event indicators.
#' @param treatment A binary vector of treatment assignments.
#' @param alternative A character value specifying the direction of the
#'   alternative hypothesis.
#'
#' @return A list containing `1 - P` for the log-rank test and an unavailable
#'   treatment-effect estimate.
#'
#' @keywords internal
#' @noRd
analyse_logrank_vectors <- function(time, event, treatment, alternative) {
  control <- treatment == 0
  lr <- logrank_test(
    groupa = time[control],
    groupb = time[!control],
    groupacensored = event[control],
    groupbcensored = event[!control]
  )
  assert_logrank_estimable(lr)

  if (alternative == "two.sided") {
    success <- 1 - lr[3]
  } else {
    # Log-rank z > 0 when control has excess events (treatment beneficial).
    # This is opposite to the Cox convention.
    # "less" => treatment beneficial => large success when z >> 0
    # "greater" => treatment harmful => large success when z << 0
    z <- lr[2]
    if (alternative == "less") {
      success <- pnorm(z)
    } else if (alternative == "greater") {
      success <- 1 - pnorm(z)
    }
  }

  list(success = success, effect = NA)
}

#' Calculate a Cox Wald result from outcome vectors
#'
#' @inheritParams analyse_logrank_vectors
#' @param h0 A finite numeric value giving the null log hazard ratio.
#'
#' @return A list containing `1 - P` for the Cox Wald test and the estimated log
#'   hazard ratio.
#'
#' @keywords internal
#' @noRd
analyse_cox_vectors <- function(time, event, treatment, alternative, h0) {
  fit_cox <- cox_wald_test_vectors_checked(time, event, treatment)
  z <- (fit_cox$estimate - h0) / fit_cox$std_error
  if (alternative == "two.sided") {
    success <- 1 - (2 * pnorm(-abs(z)))
  } else {
    # Cox z < 0 when the estimated log hazard ratio is less than h0.
    # "less" => treatment beneficial/non-inferior => large success when z << 0
    # "greater" => treatment harmful/superior to h0 => large success when z >> 0
    if (alternative == "less") {
      success <- 1 - pnorm(z)
    } else if (alternative == "greater") {
      success <- pnorm(z)
    }
  }

  list(success = success, effect = fit_cox$estimate)
}

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
analyse_bayes_surv_sufficient_stats_kernel <- function(
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

  post_lambda <- posterior_from_sufficient_stats_kernel(
    data_summ = data_summ,
    prior_surv = prior_surv,
    N_mcmc = N_mcmc,
    single_arm = single_arm,
    empty_interval = empty_interval
  )

  effect_draws <- bayes_surv_effect_draws_kernel(
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
bayes_surv_effect_draws_kernel <- function(
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

#' @title Bayesian binomial test for complete binary outcomes
#'
#' @description Updates Beta priors for one or two treatment arms and calculates
#'   the posterior probability of the specified binary-outcome hypothesis.
#'
#' @inheritParams analyse_data
#'
#' @return A list with the posterior probability of success (`success`) and
#'   posterior mean treatment effect (`effect`).
#'
#' @noRd
bayes_binomial_test <- function(
  data,
  single_arm,
  alternative,
  h0,
  prior_bin,
  bin_method,
  N_mcmc
) {
  validate_bayes_binomial_args(prior_bin, bin_method, N_mcmc)
  if (alternative == "two.sided") {
    stop(
      "Bayesian binomial analysis can only be used with alternative equal ",
      "to 'greater' or 'less'"
    )
  }

  treatment_event <- data$event[data$treatment == 1]
  if (length(treatment_event) == 0L) {
    stop("Bayesian binomial analysis requires at least one subject per arm")
  }
  control_event <- if (single_arm) {
    numeric()
  } else {
    data$event[data$treatment == 0]
  }
  if (!single_arm && length(control_event) == 0L) {
    stop("Bayesian binomial analysis requires at least one subject per arm")
  }

  bayes_binomial_count_kernel(
    events_control = sum(control_event),
    n_control = length(control_event),
    events_treatment = sum(treatment_event),
    n_treatment = length(treatment_event),
    single_arm = single_arm,
    alternative = alternative,
    h0 = h0,
    prior_bin = prior_bin,
    bin_method = bin_method,
    N_mcmc = N_mcmc
  )
}

#' Calculate a beta-binomial result from arm totals
#'
#' @description Updates the Beta prior using the number of events and subjects
#'   in each arm, then calculates the posterior probability of the specified
#'   hypothesis and the posterior mean treatment effect. These counts are the
#'   sufficient statistics for the completed binary-endpoint analysis.
#'
#' @param events_control A non-negative integer giving the number of control-arm
#'   events. Ignored in a single-arm analysis.
#' @param n_control A non-negative integer giving the number of control-arm
#'   subjects. Ignored in a single-arm analysis.
#' @param events_treatment A non-negative integer giving the number of
#'   treatment-arm events.
#' @param n_treatment A positive integer giving the number of treatment-arm
#'   subjects.
#' @inheritParams analyse_data
#'
#' @return A list containing the posterior probability of the alternative
#'   (`success`) and posterior mean event probability or treatment-control
#'   difference (`effect`). Monte Carlo analyses also retain the number of
#'   posterior draws satisfying the alternative.
#'
#' @keywords internal
#' @noRd
bayes_binomial_count_kernel <- function(
  events_control,
  n_control,
  events_treatment,
  n_treatment,
  single_arm,
  alternative,
  h0,
  prior_bin,
  bin_method,
  N_mcmc
) {
  beta_binomial_stats <- function(events, n) {
    alpha <- prior_bin[1] + events
    beta <- prior_bin[2] + n - events
    list(
      alpha = alpha,
      beta = beta,
      mean = alpha / (alpha + beta),
      variance = (alpha * beta) / ((alpha + beta)^2 * (alpha + beta + 1))
    )
  }

  beta_binomial_difference_success <- function(treatment, control) {
    integrand <- function(x) {
      treatment_density <- dbeta(x, treatment$alpha, treatment$beta)
      control_threshold <- x - h0
      if (alternative == "greater") {
        treatment_density *
          pbeta(control_threshold, control$alpha, control$beta)
      } else {
        treatment_density *
          (1 - pbeta(control_threshold, control$alpha, control$beta))
      }
    }

    integrate(integrand, lower = 0, upper = 1)$value
  }

  treatment_stats <- beta_binomial_stats(events_treatment, n_treatment)

  if (single_arm) {
    if (bin_method == "mc") {
      effect_draws <- rbeta(
        N_mcmc,
        treatment_stats$alpha,
        treatment_stats$beta
      )
      success_indicator <- if (alternative == "greater") {
        effect_draws > h0
      } else {
        effect_draws < h0
      }
      success <- mean(success_indicator)
      effect <- mean(effect_draws)
    } else {
      effect <- treatment_stats$mean
      if (bin_method == "normal") {
        effect_se <- sqrt(treatment_stats$variance)
        success <- if (alternative == "greater") {
          1 - pnorm(h0, mean = effect, sd = effect_se)
        } else {
          pnorm(h0, mean = effect, sd = effect_se)
        }
      } else {
        success <- if (alternative == "greater") {
          1 - pbeta(h0, treatment_stats$alpha, treatment_stats$beta)
        } else {
          pbeta(h0, treatment_stats$alpha, treatment_stats$beta)
        }
      }
    }
    result <- list(success = success, effect = effect)
    if (bin_method == "mc") {
      result <- set_analysis_mc_counts(
        result,
        successes = sum(success_indicator),
        draws = length(success_indicator)
      )
    }
    return(result)
  }

  control_stats <- beta_binomial_stats(events_control, n_control)

  if (bin_method == "mc") {
    effect_draws <-
      rbeta(N_mcmc, treatment_stats$alpha, treatment_stats$beta) -
      rbeta(N_mcmc, control_stats$alpha, control_stats$beta)
    success_indicator <- if (alternative == "greater") {
      effect_draws > h0
    } else {
      effect_draws < h0
    }
    success <- mean(success_indicator)
    effect <- mean(effect_draws)
  } else if (bin_method == "normal") {
    effect <- treatment_stats$mean - control_stats$mean
    effect_se <- sqrt(treatment_stats$variance + control_stats$variance)
    if (alternative == "greater") {
      success <- 1 - pnorm(h0, mean = effect, sd = effect_se)
    } else if (alternative == "less") {
      success <- pnorm(h0, mean = effect, sd = effect_se)
    }
  } else if (bin_method == "quadrature") {
    effect <- treatment_stats$mean - control_stats$mean
    success <- beta_binomial_difference_success(treatment_stats, control_stats)
  }

  result <- list(success = success, effect = effect)
  if (bin_method == "mc") {
    result <- set_analysis_mc_counts(
      result,
      successes = sum(success_indicator),
      draws = length(success_indicator)
    )
  }
  result
}

#' @title Validate complete binary outcomes
#'
#' @description Ensures a binary endpoint analysis is supplied only event
#'   indicators and that censored subjects have complete follow-up.
#'
#' @noRd
assert_complete_binary_outcomes <- function(data, end_of_study, method_label) {
  tryCatch(
    validate_binary_vector(data$event, "data$event"),
    error = function(error) {
      stop(
        method_label,
        " analysis requires binary event outcomes: ",
        conditionMessage(error)
      )
    }
  )
  validate_nonnegative_numeric_vector(data$time, "data$time")
  if (any(data$event == 0 & data$time < end_of_study)) {
    stop(
      method_label,
      " analysis requires all censored subjects to be followed ",
      "to 'end_of_study' or imputed before analysis"
    )
  }
}

#' @title Estimate a frequentist binary risk difference
#'
#' @description Calculates the treatment-minus-control difference in event
#'   proportions and its unpooled binomial variance.
#'
#' @inheritParams analyse_data
#'
#' @return A list containing the risk-difference estimate and variance.
#'
#' @noRd
risk_difference_estimate_checked <- function(data, end_of_study) {
  assert_complete_binary_outcomes(data, end_of_study, "Risk-difference")

  treatment_event <- data$event[data$treatment == 1]
  control_event <- data$event[data$treatment == 0]
  if (length(treatment_event) == 0 || length(control_event) == 0) {
    stop("Risk-difference analysis requires at least one subject per arm")
  }

  risk_difference_estimate_count_kernel(
    events_control = sum(control_event),
    n_control = length(control_event),
    events_treatment = sum(treatment_event),
    n_treatment = length(treatment_event)
  )
}

#' Estimate a binary risk difference from arm totals
#'
#' @param events_control A non-negative integer giving the number of control-arm
#'   events.
#' @param n_control A positive integer giving the number of control-arm
#'   subjects.
#' @param events_treatment A non-negative integer giving the number of
#'   treatment-arm events.
#' @param n_treatment A positive integer giving the number of treatment-arm
#'   subjects.
#'
#' @return A list containing the treatment-minus-control event-risk difference
#'   (`estimate`) and its unpooled binomial variance (`variance`).
#'
#' @keywords internal
#' @noRd
risk_difference_estimate_count_kernel <- function(
  events_control,
  n_control,
  events_treatment,
  n_treatment
) {
  treatment_risk <- events_treatment / n_treatment
  control_risk <- events_control / n_control
  estimate <- treatment_risk - control_risk
  variance <- treatment_risk *
    (1 - treatment_risk) /
    n_treatment +
    control_risk * (1 - control_risk) / n_control

  if (!is.finite(estimate) || !is.finite(variance) || variance < 0) {
    stop(
      "Risk-difference analysis is non-estimable: the estimate or variance is ",
      "not finite and non-negative."
    )
  }

  list(estimate = estimate, variance = variance)
}

#' @title Test a frequentist binary risk difference
#'
#' @description Performs a Wald test of the treatment-minus-control event-risk
#'   difference against the supplied null value.
#'
#' @inheritParams analyse_data
#'
#' @return A list containing the success score, risk-difference estimate,
#'   variance, and standard error.
#'
#' @noRd
risk_difference_wald_test_checked <- function(
  data,
  end_of_study,
  alternative,
  h0
) {
  fit <- risk_difference_estimate_checked(data, end_of_study)
  risk_difference_wald_from_estimate_kernel(fit, alternative, h0)
}

#' Test a binary risk difference from arm totals
#'
#' @inheritParams risk_difference_estimate_count_kernel
#' @param alternative A character value giving the direction of the alternative
#'   hypothesis: `"less"`, `"greater"`, or `"two.sided"`.
#' @param h0 A finite numeric value giving the null treatment-minus-control
#'   event-risk difference.
#'
#' @return A list containing the success score (`success`), estimated risk
#'   difference (`estimate`), variance, and standard error.
#'
#' @keywords internal
#' @noRd
risk_difference_wald_count_kernel <- function(
  events_control,
  n_control,
  events_treatment,
  n_treatment,
  alternative,
  h0
) {
  fit <- risk_difference_estimate_count_kernel(
    events_control = events_control,
    n_control = n_control,
    events_treatment = events_treatment,
    n_treatment = n_treatment
  )
  risk_difference_wald_from_estimate_kernel(fit, alternative, h0)
}

#' Test a binary risk difference from its estimate and variance
#'
#' @param fit A list containing a finite risk-difference estimate and a
#'   non-negative variance.
#' @inheritParams risk_difference_wald_count_kernel
#'
#' @return A list containing the success score (`success`), estimated risk
#'   difference (`estimate`), variance, and standard error.
#'
#' @keywords internal
#' @noRd
risk_difference_wald_from_estimate_kernel <- function(fit, alternative, h0) {
  if (fit$variance == 0) {
    stop(
      "Risk-difference analysis is non-estimable: the estimated variance is ",
      "zero because both arms have no outcome variation."
    )
  }

  std_error <- sqrt(fit$variance)
  statistic <- (fit$estimate - h0) / std_error
  success <- switch(
    alternative,
    "less" = 1 - pnorm(statistic),
    "greater" = pnorm(statistic),
    "two.sided" = 1 - 2 * pnorm(-abs(statistic))
  )

  list(
    success = success,
    estimate = fit$estimate,
    variance = fit$variance,
    std_error = std_error
  )
}

#' @title Calculate a two-sample log-rank test
#'
#' @description Calculates a log-rank test from the follow-up times and event
#'   indicators in two independent treatment groups.
#'
#' @param groupa A numeric vector of non-negative follow-up times for group A.
#' @param groupb A numeric vector of non-negative follow-up times for group B.
#' @param groupacensored A zero-one integer vector indicating events in group A,
#'   where `1` denotes an event and `0` denotes censoring.
#' @param groupbcensored A zero-one integer vector indicating events in group B,
#'   where `1` denotes an event and `0` denotes censoring.
#' @param onlyz A single logical value indicating whether to return only the
#'   standardized log-rank statistic. The default is `FALSE`.
#'
#' @return A numeric vector containing the chi-squared statistic, standardized
#'   statistic, and two-sided *P*-value. When `onlyz = TRUE`, only the
#'   standardized statistic is returned.
#'
#' @examples
#' T1 <- c(6, 6, 6, 6, 7, 9, 10, 10, 11, 13, 16, 17, 19, 20, 22, 23, 25, 32, 32, 34, 35)
#' E1 <- c(1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0)
#' T2 <- c(1, 1, 2, 2, 3, 4, 4, 5, 5, 8, 8, 8, 8, 11, 11, 12, 12, 15, 17, 22, 23)
#' E2 <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' logrank_test(T1, T2, E1, E2)
#' #1.679294e+01 -4.097919e+00, 4.168809e-05
#'
#' @noRd
#' @keywords internal
logrank_test <- function(
  groupa,
  groupb,
  groupacensored,
  groupbcensored,
  onlyz = FALSE
) {
  logrank_instance(groupa, groupb, groupacensored, groupbcensored, onlyz)
}

#' @title Assert that a log-rank result is estimable
#'
#' @description Stops with a diagnostic error when a log-rank statistic cannot
#'   be computed as a finite comparison between treatment groups.
#'
#' @noRd
assert_logrank_estimable <- function(lr) {
  if (
    length(lr) < 3 ||
      any(is.na(lr[1:3])) ||
      any(!is.finite(lr[1:3]))
  ) {
    stop(
      "Log-rank analysis is non-estimable: the test statistic is not finite. ",
      "This can occur when there are no events or insufficient information ",
      "to compare treatment groups."
    )
  }

  invisible(TRUE)
}

#' @title Calculate a Cox proportional hazards Wald test
#'
#' @description Fits a Cox proportional hazards model for the treatment effect
#'   and converts fitting warnings into clear non-estimability errors.
#'
#' @inheritParams analyse_data
#'
#' @return A list with the estimated log hazard ratio (`estimate`) and its
#'   standard error (`std_error`).
#'
#' @noRd
cox_wald_test_checked <- function(data, ...) {
  cox_wald_test_vectors_checked(
    time = data$time,
    event = data$event,
    treatment = data$treatment,
    ...
  )
}

#' Calculate a checked Cox Wald test from outcome vectors
#'
#' @inheritParams analyse_logrank_vectors
#'
#' @return A list with the estimated log hazard ratio (`estimate`) and its
#'   standard error (`std_error`).
#'
#' @keywords internal
#' @noRd
cox_wald_test_vectors_checked <- function(time, event, treatment, ...) {
  fit_state <- new.env(parent = emptyenv())
  fit_state$warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(
      cox_wald_test_vectors(time, event, treatment, ...),
      warning = function(w) {
        fit_state$warnings <- c(
          fit_state$warnings,
          conditionMessage(w)
        )
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      stop(
        "Cox analysis is non-estimable: the Cox fitter failed: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  if (length(fit_state$warnings) > 0L) {
    stop(
      "Cox analysis is non-estimable: the Cox model did not produce a ",
      "reliable Wald statistic. The Cox fitter reported: ",
      paste(unique(fit_state$warnings), collapse = " | "),
      call. = FALSE
    )
  }

  assert_cox_estimable(fit)
  fit
}

#' @title Assert that a Cox result is estimable
#'
#' @description Ensures a Cox treatment-effect estimate and its standard error
#'   are finite before they are used in an adaptive decision.
#'
#' @noRd
assert_cox_estimable <- function(fit) {
  if (
    length(fit$estimate) != 1 ||
      length(fit$std_error) != 1 ||
      is.na(fit$estimate) ||
      is.na(fit$std_error) ||
      !is.finite(fit$estimate) ||
      !is.finite(fit$std_error) ||
      fit$std_error <= 0
  ) {
    stop(
      "Cox analysis is non-estimable: the treatment effect or standard ",
      "error is not finite and positive. This can occur when there are no ",
      "events, separation, or insufficient information to estimate the ",
      "treatment effect."
    )
  }

  invisible(TRUE)
}

#' @title Fit a Cox Wald test using an available survival calculation
#'
#' @description Uses the lower-overhead `coxph.fit()` calculation when its
#'   arguments are compatible with the installed version of `survival`, and
#'   otherwise uses [survival::coxph()]. The `engine` and `compatibility`
#'   arguments support equivalence testing of the two calculations.
#'
#' @noRd
cox_wald_test <- function(
  data,
  engine = c("auto", "fast", "public"),
  compatibility = coxph_fit_compatibility()
) {
  cox_wald_test_vectors(
    time = data$time,
    event = data$event,
    treatment = data$treatment,
    engine = engine,
    compatibility = compatibility
  )
}

#' Fit a Cox Wald test from outcome vectors
#'
#' @inheritParams analyse_logrank_vectors
#' @param engine A character value selecting the automatic, faster, or exported
#'   Cox calculation.
#' @param compatibility The result of `coxph_fit_compatibility()`.
#'
#' @return A list with the estimated log hazard ratio (`estimate`) and its
#'   standard error (`std_error`).
#'
#' @keywords internal
#' @noRd
cox_wald_test_vectors <- function(
  time,
  event,
  treatment,
  engine = c("auto", "fast", "public"),
  compatibility = coxph_fit_compatibility()
) {
  engine <- match.arg(engine)
  if (engine == "public") {
    return(cox_wald_test_public_vectors(time, event, treatment))
  }

  fast_available <- isTRUE(compatibility$compatible) &&
    is.function(compatibility$fitter)
  if (!fast_available) {
    if (engine == "fast") {
      stop(
        "The survival::coxph.fit() fast path is unavailable: ",
        compatibility$reason,
        call. = FALSE
      )
    }
    return(cox_wald_test_public_vectors(time, event, treatment))
  }

  fast_fit <- tryCatch(
    cox_wald_test_fast_vectors(time, event, treatment, compatibility$fitter),
    error = identity
  )
  if (!inherits(fast_fit, "error")) {
    return(fast_fit)
  }
  if (engine == "fast") {
    stop(fast_fit)
  }

  cox_wald_test_public_vectors(time, event, treatment)
}

#' @title Check availability of the lower-overhead Cox calculation
#'
#' @description Records the installed `survival` version and verifies every
#'   unexported fitter argument used by the package. Installations outside the
#'   supported compatibility boundary use the exported fallback.
#'
#' @noRd
coxph_fit_compatibility <- local({
  cached <- NULL

  function(refresh = FALSE) {
    if (!refresh && !is.null(cached)) {
      return(cached)
    }

    version <- as.character(utils::packageVersion("survival"))
    fitter <- tryCatch(
      getFromNamespace("coxph.fit", "survival"),
      error = function(e) NULL
    )
    required_arguments <- c(
      "x",
      "y",
      "strata",
      "offset",
      "init",
      "control",
      "weights",
      "method",
      "rownames",
      "resid",
      "nocenter"
    )
    version_supported <- utils::compareVersion(version, "3.2-0") >= 0L
    missing_arguments <- if (is.function(fitter)) {
      setdiff(required_arguments, names(formals(fitter)))
    } else {
      required_arguments
    }
    compatible <- version_supported &&
      is.function(fitter) &&
      length(missing_arguments) == 0L
    reason <- if (!version_supported) {
      paste0("survival ", version, " predates the guarded fast path")
    } else if (!is.function(fitter)) {
      "survival::coxph.fit() could not be resolved"
    } else if (length(missing_arguments) > 0L) {
      paste0(
        "survival::coxph.fit() lacks required argument(s): ",
        paste(missing_arguments, collapse = ", ")
      )
    } else {
      paste0("compatible survival ", version, " signature")
    }

    cached <<- list(
      compatible = compatible,
      fitter = if (compatible) fitter else NULL,
      version = version,
      reason = reason
    )
    cached
  }
})

#' Fit a Cox Wald test from outcome vectors with `coxph.fit()`
#'
#' @inheritParams analyse_logrank_vectors
#' @param fitter A compatible `survival::coxph.fit()` function.
#'
#' @return A list with the estimated log hazard ratio (`estimate`) and its
#'   standard error (`std_error`).
#'
#' @keywords internal
#' @noRd
cox_wald_test_fast_vectors <- function(time, event, treatment, fitter) {
  y <- survival::Surv(time, event)
  x <- matrix(
    as.double(treatment),
    ncol = 1L,
    dimnames = list(NULL, "treatment")
  )
  fit <- fitter(
    x = x,
    y = y,
    strata = NULL,
    offset = NULL,
    init = NULL,
    control = survival::coxph.control(),
    weights = NULL,
    method = "efron",
    rownames = NULL,
    resid = FALSE,
    nocenter = NULL
  )

  valid_result <- is.list(fit) &&
    is.numeric(fit$coefficients) &&
    identical(names(fit$coefficients), "treatment") &&
    is.matrix(fit$var) &&
    identical(dim(fit$var), c(1L, 1L))
  if (!valid_result) {
    stop(
      "survival::coxph.fit() returned an incompatible result structure",
      call. = FALSE
    )
  }

  list(
    estimate = unname(fit$coefficients[["treatment"]]),
    std_error = sqrt(unname(fit$var[1L, 1L]))
  )
}

#' Fit a Cox Wald test from outcome vectors with [survival::coxph()]
#'
#' @inheritParams analyse_logrank_vectors
#'
#' @return A list with the estimated log hazard ratio (`estimate`) and its
#'   standard error (`std_error`).
#'
#' @keywords internal
#' @noRd
cox_wald_test_public_vectors <- function(time, event, treatment) {
  data <- data.frame(
    time = time,
    event = event,
    treatment = treatment
  )
  fit <- survival::coxph(
    survival::Surv(time, event) ~ treatment,
    data = data,
    ties = "efron",
    singular.ok = FALSE,
    model = FALSE,
    x = FALSE,
    y = FALSE
  )
  estimate <- stats::coef(fit)["treatment"]
  variance <- stats::vcov(fit)["treatment", "treatment"]

  list(
    estimate = unname(estimate),
    std_error = sqrt(unname(variance))
  )
}
