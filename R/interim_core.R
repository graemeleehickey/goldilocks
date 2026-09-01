#' Evaluate one prepared interim analysis
#'
#' @description Runs the posterior-predictive calculation shared by simulated
#'   trials and externally observed interim data. Callers are responsible for
#'   preparing the internal enrollment and imputation indicators.
#'
#' @param data_interim Data frame containing `time`, `event`, `treatment`,
#'   `subject_enrolled`, `subject_impute_success`, and
#'   `subject_impute_futility`.
#' @param look Integer index of the interim look.
#' @param planned_N Planned enrollment at the look.
#' @param calendar_time Calendar time of the data cut.
#' @param active_followup Number of enrolled subjects still under follow-up.
#' @param check_futility Whether to calculate success at the maximum sample
#'   size.
#' @inheritParams survival_adapt
#'
#' @return A list containing predictive probabilities, the decision, Monte
#'   Carlo summaries, diagnostics, and a one-row decision trace.
#'
#' @keywords internal
#' @noRd
evaluate_interim_core <- function(
  data_interim,
  look,
  planned_N,
  calendar_time,
  active_followup,
  end_of_study,
  cutpoints,
  single_arm,
  prior_surv,
  prior_bin,
  bin_method,
  alternative,
  h0,
  Fn,
  Sn,
  prob_ha,
  N_impute,
  N_mcmc,
  mc_conf_level,
  empty_interval,
  method,
  binary_imputation,
  check_futility
) {
  required_columns <- c(
    "time",
    "event",
    "treatment",
    "subject_enrolled",
    "subject_impute_success",
    "subject_impute_futility"
  )
  missing_columns <- setdiff(required_columns, names(data_interim))
  if (!is.data.frame(data_interim) || length(missing_columns) > 0L) {
    stop(
      "'data_interim' must be a data frame containing columns: ",
      paste(required_columns, collapse = ", ")
    )
  }

  data <- data_interim[
    data_interim$subject_enrolled,
    c("time", "event", "treatment"),
    drop = FALSE
  ]

  warning_state <- new.env(parent = emptyenv())
  warning_state$messages <- character()
  warning_state$empty_interval_fallbacks <- character()
  capture_warning <- function(warning) {
    warning_state$messages <- unique(c(
      warning_state$messages,
      conditionMessage(warning)
    ))
  }
  capture_empty_interval <- function(condition) {
    detail <- condition$details
    entries <- paste0(
      condition$policy,
      ": treatment=",
      detail$treatment,
      ", interval=",
      detail$interval
    )
    warning_state$empty_interval_fallbacks <- unique(c(
      warning_state$empty_interval_fallbacks,
      entries
    ))
  }

  data_summ <- posterior_sufficient_stats(
    data = data,
    cutpoints = cutpoints,
    single_arm = single_arm
  )
  post_lambda <- withCallingHandlers(
    posterior_from_sufficient_stats(
      data_summ = data_summ,
      prior_surv = prior_surv,
      N_mcmc = N_impute,
      single_arm = single_arm,
      empty_interval = empty_interval
    ),
    warning = capture_warning,
    goldilocks_empty_interval = capture_empty_interval
  )
  posterior_diagnostics <- gamma_posterior_diagnostics(
    data_summ = data_summ,
    prior_surv = prior_surv,
    cutpoints = cutpoints,
    end_of_study = end_of_study,
    single_arm = single_arm,
    empty_interval = empty_interval
  )
  interval_widths <- if (method == "bayes-surv") {
    endpoint_interval_widths(cutpoints, end_of_study)
  } else {
    NULL
  }

  maximum_successes <- 0L
  current_successes <- 0L
  inner_mc_uncertain_now <- 0L
  inner_mc_uncertain_max <- 0L
  for (j in seq_len(N_impute)) {
    hazard <- post_lambda[j, , , drop = FALSE]

    stop_check <- withCallingHandlers(
      test_stop_success(
        data = data_interim,
        hazard = hazard,
        end_of_study = end_of_study,
        cutpoints = cutpoints,
        interval_widths = interval_widths,
        single_arm = single_arm,
        prior_surv = prior_surv,
        N_mcmc = N_mcmc,
        method = method,
        alternative = alternative,
        h0 = h0,
        prior_bin = prior_bin,
        bin_method = bin_method,
        binary_imputation = if (method %in% c("bayes-bin", "riskdiff")) {
          binary_imputation
        } else {
          "event-time"
        },
        empty_interval = empty_interval,
        check_futility = check_futility,
        mc_conf_level = mc_conf_level,
        prob_ha = prob_ha
      ),
      warning = capture_warning,
      goldilocks_empty_interval = capture_empty_interval
    )

    analysis_now <- stop_check$classification_now
    current_successes <- current_successes + as.integer(analysis_now$crossed)
    inner_mc_uncertain_now <- inner_mc_uncertain_now +
      as.integer(analysis_now$uncertain)

    if (check_futility) {
      analysis_max <- stop_check$classification_max
      maximum_successes <- maximum_successes +
        as.integer(analysis_max$crossed)
      inner_mc_uncertain_max <- inner_mc_uncertain_max +
        as.integer(analysis_max$uncertain)
    }
  }

  ppp_now_mc <- monte_carlo_probability_summary(
    successes = current_successes,
    draws = N_impute,
    threshold = Sn,
    direction = "greater",
    confidence = mc_conf_level
  )
  ppp_max_mc <- if (check_futility) {
    monte_carlo_probability_summary(
      successes = maximum_successes,
      draws = N_impute,
      threshold = Fn,
      direction = "less",
      confidence = mc_conf_level
    )
  } else {
    NULL
  }

  decision <- if (ppp_now_mc$point_crossed) {
    "stop_expected_success"
  } else if (check_futility && ppp_max_mc$point_crossed) {
    "stop_futility"
  } else {
    "continue"
  }
  decision_reason <- if (decision == "stop_expected_success") {
    "expected_success_estimate_above_threshold"
  } else if (decision == "stop_futility") {
    "futility_estimate_below_threshold"
  } else {
    "continue_thresholds_not_crossed"
  }

  trace <- data.frame(
    look = as.integer(look),
    planned_N = as.integer(planned_N),
    calendar_time = calendar_time,
    active_followup = as.integer(active_followup),
    N_enrolled = nrow(data),
    N_treatment = sum(data$treatment == 1),
    N_control = sum(data$treatment == 0),
    events_treatment = sum(data$event[data$treatment == 1]),
    events_control = sum(data$event[data$treatment == 0]),
    N_pending = sum(
      data_interim$subject_enrolled &
        data_interim$subject_impute_success
    ),
    N_not_enrolled = sum(data_interim$subject_impute_futility),
    ppp_stop_now = ppp_now_mc$estimate,
    ppp_stop_now_mcse = ppp_now_mc$mcse,
    ppp_stop_now_lower = ppp_now_mc$lower,
    ppp_stop_now_upper = ppp_now_mc$upper,
    ppp_stop_now_draws = ppp_now_mc$draws,
    success_threshold = Sn,
    ppp_success_at_max = if (check_futility) {
      ppp_max_mc$estimate
    } else {
      NA_real_
    },
    ppp_success_at_max_mcse = if (check_futility) {
      ppp_max_mc$mcse
    } else {
      NA_real_
    },
    ppp_success_at_max_lower = if (check_futility) {
      ppp_max_mc$lower
    } else {
      NA_real_
    },
    ppp_success_at_max_upper = if (check_futility) {
      ppp_max_mc$upper
    } else {
      NA_real_
    },
    ppp_success_at_max_draws = if (check_futility) {
      ppp_max_mc$draws
    } else {
      NA_integer_
    },
    futility_threshold = if (check_futility) Fn else NA_real_,
    inner_mc_uncertain_stop_now = inner_mc_uncertain_now,
    inner_mc_uncertain_success_at_max = if (check_futility) {
      inner_mc_uncertain_max
    } else {
      NA_integer_
    },
    decision = decision,
    decision_reason = decision_reason,
    empty_interval_fallback_count = length(
      warning_state$empty_interval_fallbacks
    ),
    empty_interval_fallbacks = paste(
      warning_state$empty_interval_fallbacks,
      collapse = " | "
    ),
    warning_count = length(warning_state$messages),
    warning_messages = paste(warning_state$messages, collapse = " | "),
    stringsAsFactors = FALSE
  )

  monte_carlo <- data.frame(
    estimand = c(
      "success_if_stop_now",
      if (check_futility) "success_at_maximum" else character()
    ),
    successes = c(
      ppp_now_mc$successes,
      if (check_futility) ppp_max_mc$successes else integer()
    ),
    draws = c(
      ppp_now_mc$draws,
      if (check_futility) ppp_max_mc$draws else integer()
    ),
    estimate = c(
      ppp_now_mc$estimate,
      if (check_futility) ppp_max_mc$estimate else numeric()
    ),
    mcse = c(
      ppp_now_mc$mcse,
      if (check_futility) ppp_max_mc$mcse else numeric()
    ),
    lower = c(
      ppp_now_mc$lower,
      if (check_futility) ppp_max_mc$lower else numeric()
    ),
    upper = c(
      ppp_now_mc$upper,
      if (check_futility) ppp_max_mc$upper else numeric()
    ),
    threshold = c(
      ppp_now_mc$threshold,
      if (check_futility) ppp_max_mc$threshold else numeric()
    ),
    direction = c(
      ppp_now_mc$direction,
      if (check_futility) ppp_max_mc$direction else character()
    ),
    point_crossed = c(
      ppp_now_mc$point_crossed,
      if (check_futility) ppp_max_mc$point_crossed else logical()
    ),
    bound_crossed = c(
      ppp_now_mc$bound_crossed,
      if (check_futility) ppp_max_mc$bound_crossed else logical()
    ),
    stringsAsFactors = FALSE
  )

  list(
    ppp_success = ppp_now_mc$estimate,
    ppp_success_at_max = if (check_futility) {
      ppp_max_mc$estimate
    } else {
      NA_real_
    },
    decision = decision,
    decision_reason = decision_reason,
    monte_carlo = monte_carlo,
    diagnostics = list(
      inner_mc_uncertain_stop_now = inner_mc_uncertain_now,
      inner_mc_uncertain_success_at_max = if (check_futility) {
        inner_mc_uncertain_max
      } else {
        NA_integer_
      },
      empty_interval_fallbacks = warning_state$empty_interval_fallbacks,
      warnings = warning_state$messages,
      prior = gamma_prior_diagnostics(
        prior_surv = prior_surv,
        cutpoints = cutpoints,
        end_of_study = end_of_study,
        single_arm = single_arm,
        stage = "interim"
      ),
      posterior = posterior_diagnostics
    ),
    trace = trace
  )
}

#' Bind interim posterior diagnostic rows
#'
#' @param rows List containing one posterior diagnostic data frame per look.
#'
#' @return A stable data frame with one row per look, arm, and interval.
#'
#' @keywords internal
#' @noRd
new_interim_posterior_diagnostics <- function(rows) {
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) > 0L) {
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    return(out)
  }

  data.frame(
    look = integer(),
    planned_N = integer(),
    calendar_time = numeric(),
    arm = character(),
    interval = integer(),
    interval_start = numeric(),
    interval_end = numeric(),
    exposed_subjects = integer(),
    observed_exposure = numeric(),
    observed_events = integer(),
    empty_interval = logical(),
    empty_interval_policy = character(),
    effective_exposure = numeric(),
    effective_events = numeric(),
    prior_shape = numeric(),
    prior_rate = numeric(),
    posterior_shape = numeric(),
    posterior_rate = numeric(),
    posterior_mean_hazard = numeric(),
    posterior_sd_hazard = numeric(),
    stringsAsFactors = FALSE
  )
}
