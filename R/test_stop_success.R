#' @title Evaluate success for one predictively completed trial
#'
#' @description Constructs one completed trial from the predictive imputations
#'   and evaluates the prespecified success criterion at the current and, when
#'   requested, maximum sample size.
#'
#' @inheritParams survival_adapt
#' @inheritParams sim_comp_data
#' @param imputations A list of interim predictive imputations returned by
#'   `impute_predictive_draws()`.
#' @param draw A positive integer identifying the predictive imputation to
#'   analyze.
#' @param check_futility A single logical value indicating whether the design
#'   includes a futility assessment at the maximum sample size.
#' @param interval_widths A numeric vector giving the analysis-interval
#'   durations through `end_of_study`, used for the Bayesian survival analysis.
#' @param prob_ha A single numeric probability in `[0, 1]` defining success in
#'   the completed-data analysis.
#' @param mc_conf_level A single numeric probability strictly between `0.5` and
#'   `1`, giving the level of diagnostic finite Monte Carlo bounds.
#'
#' @return A numeric vector containing indicators of success at the current and
#'   maximum sample sizes and indicators of material inner Monte Carlo
#'   uncertainty.
#' @noRd
test_stop_success <- function(
  data,
  imputations,
  draw,
  end_of_study,
  cutpoints,
  interval_widths,
  single_arm,
  prior_surv,
  N_mcmc,
  method,
  alternative,
  h0,
  prior_bin,
  bin_method,
  empty_interval,
  check_futility,
  prob_ha,
  mc_conf_level
) {
  ##############################################################################
  ### Test for success at current sample size (-> stop for success)
  ##############################################################################

  data_success_impute <- materialize_predictive_draw(
    data_in = data,
    imputations = imputations,
    draw = draw
  )

  if (method == "bayes-surv") {
    # This is stage two of the posterior-predictive calculation. `hazard`
    # above came from the observed-data posterior and generated one completed
    # trial. We now summarize that completed trial and form a *fresh* posterior
    # from its statistics plus the original prior.
    #
    # `rows` restricts the current-sample-size analysis to enrolled subjects
    # without first copying time/event/treatment into another data frame.
    data_summ_now <- posterior_sufficient_stats(
      data = data_success_impute,
      cutpoints = cutpoints,
      single_arm = single_arm,
      rows = data_success_impute$subject_enrolled
    )
    success_now <- analyse_bayes_surv_sufficient_stats_kernel(
      data_summ = data_summ_now,
      cutpoints = cutpoints,
      end_of_study = end_of_study,
      prior_surv = prior_surv,
      N_mcmc = N_mcmc,
      single_arm = single_arm,
      alternative = alternative,
      h0 = h0,
      empty_interval = empty_interval,
      interval_widths = interval_widths
    )
  } else {
    # Frequentist and binary-Bayesian methods still require patient-level
    # outcomes, so retain their existing analysis data frame.
    time <- NULL
    event <- NULL
    subject_enrolled <- NULL
    data_now <- subset(
      data_success_impute,
      subset = subject_enrolled,
      select = c(time, event, treatment)
    )
    success_now <- analyse_data(
      data = data_now,
      cutpoints = cutpoints,
      end_of_study = end_of_study,
      prior_surv = prior_surv,
      N_mcmc = N_mcmc,
      single_arm = single_arm,
      method = method,
      alternative = alternative,
      h0 = h0,
      prior_bin = prior_bin,
      bin_method = bin_method,
      empty_interval = empty_interval
    )
  }

  ##############################################################################
  ### Test for success at maximum sample size (-> stop for futility)
  ##############################################################################

  if (check_futility) {
    # Take the already imputed data for expected success and add the batched
    # outcomes for subjects not yet enrolled.
    data_futility_impute <- materialize_predictive_draw(
      data_in = data_success_impute,
      imputations = imputations,
      draw = draw,
      include_future = TRUE
    )

    if (method == "bayes-surv") {
      # The maximum-sample-size data include both the success imputations for
      # enrolled subjects and the futility imputations for future subjects.
      # Summarize all rows because this branch represents a completed maximum
      # trial rather than the currently enrolled cohort.
      data_summ_max <- posterior_sufficient_stats(
        data = data_futility_impute,
        cutpoints = cutpoints,
        single_arm = single_arm
      )
      success_max <- analyse_bayes_surv_sufficient_stats_kernel(
        data_summ = data_summ_max,
        cutpoints = cutpoints,
        end_of_study = end_of_study,
        prior_surv = prior_surv,
        N_mcmc = N_mcmc,
        single_arm = single_arm,
        alternative = alternative,
        h0 = h0,
        empty_interval = empty_interval,
        interval_widths = interval_widths
      )
    } else {
      # Create the patient-level data required by the other analysis methods.
      time <- NULL
      event <- NULL
      treatment <- NULL
      data_max <- subset(
        data_futility_impute,
        select = c(time, event, treatment)
      )
      success_max <- analyse_data(
        data = data_max,
        cutpoints = cutpoints,
        end_of_study = end_of_study,
        prior_surv = prior_surv,
        N_mcmc = N_mcmc,
        single_arm = single_arm,
        method = method,
        alternative = alternative,
        h0 = h0,
        prior_bin = prior_bin,
        bin_method = bin_method,
        empty_interval = empty_interval
      )
    }
  } else {
    success_max <- NA
  }

  classification_now <- classify_completed_analysis(
    success_now,
    prob_ha = prob_ha,
    mc_conf_level = mc_conf_level
  )
  classification_max <- if (check_futility) {
    classify_completed_analysis(
      success_max,
      prob_ha = prob_ha,
      mc_conf_level = mc_conf_level
    )
  } else {
    NULL
  }

  return(list(
    success_now = success_now,
    success_max = success_max,
    classification_now = classification_now,
    classification_max = classification_max
  ))
}

#' Apply the binary success rule across predictive replicates
#'
#' @description Applies the prespecified risk-difference or beta-binomial
#'   analysis to the arm-specific event counts for every predictive replicate.
#'   Identical sufficient statistics have identical results for deterministic
#'   analyses and are therefore evaluated once within a look. A Monte Carlo
#'   beta-binomial analysis uses fresh posterior draws for every replicate,
#'   including replicates with identical event counts.
#'
#' @param count_states A list returned by
#'   `predictive_binary_count_states()`, containing arm-specific event counts
#'   and sample sizes for the current and maximum-sample completions.
#' @param single_arm A single logical value indicating a one-arm design.
#' @param N_mcmc A positive integer giving the number of beta-posterior draws
#'   for `method = "bayes-bin"` and `bin_method = "mc"`.
#' @param method A character value specifying `"riskdiff"` or `"bayes-bin"`.
#' @param alternative A character value giving the direction of the alternative
#'   hypothesis.
#' @param h0 A finite numeric value giving the null event-probability difference
#'   or, in a single-arm analysis, the null event probability.
#' @param prior_bin A length-two positive numeric vector giving the Beta prior
#'   shape parameters for events and non-events.
#' @param bin_method A character value specifying `"mc"`, `"normal"`, or
#'   `"quadrature"` for a Bayesian binary analysis.
#' @param check_futility A single logical value indicating whether to assess
#'   success after continuing to the maximum sample size.
#' @param prob_ha A numeric probability in `[0, 1]` defining success for a
#'   completed replicate.
#' @param mc_conf_level A numeric probability strictly between `0.5` and `1`
#'   giving the level of diagnostic Monte Carlo bounds.
#'
#' @return A list containing the numbers of successful current and maximum-
#'   sample predictive replicates, counts of classifications with material
#'   inner Monte Carlo uncertainty, and a summary of repeated sufficient-
#'   statistic combinations.
#'
#' @keywords internal
#' @noRd
analyse_predictive_binary_counts <- function(
  count_states,
  single_arm,
  N_mcmc,
  method,
  alternative,
  h0,
  prior_bin,
  bin_method,
  check_futility,
  prob_ha,
  mc_conf_level
) {
  n_draws <- nrow(count_states$current)
  if (check_futility) {
    combined <- rbind(count_states$current, count_states$maximum)
    draw_order <- as.vector(rbind(
      seq_len(n_draws),
      n_draws + seq_len(n_draws)
    ))
    combined <- combined[draw_order, , drop = FALSE]
    stage <- rep(c("current", "maximum"), n_draws)
  } else {
    combined <- count_states$current
    stage <- rep.int("current", n_draws)
  }
  rownames(combined) <- NULL

  keys <- binary_count_analysis_keys(
    states = combined,
    method = method,
    single_arm = single_arm,
    alternative = alternative,
    h0 = h0,
    prior_bin = prior_bin,
    bin_method = bin_method,
    N_mcmc = N_mcmc
  )
  deterministic <- method == "riskdiff" || bin_method != "mc"
  unique_state <- !duplicated(keys)
  analysis_rows <- if (deterministic) {
    which(unique_state)
  } else {
    seq_len(nrow(combined))
  }

  analyses_unique <- lapply(analysis_rows, function(row) {
    analyse_binary_count_state_kernel(
      state = combined[row, , drop = FALSE],
      single_arm = single_arm,
      N_mcmc = N_mcmc,
      method = method,
      alternative = alternative,
      h0 = h0,
      prior_bin = prior_bin,
      bin_method = bin_method
    )
  })
  analyses <- if (deterministic) {
    analyses_unique[match(keys, keys[unique_state])]
  } else {
    analyses_unique
  }
  classifications <- lapply(analyses, function(analysis) {
    classify_completed_analysis(
      analysis,
      prob_ha = prob_ha,
      mc_conf_level = mc_conf_level
    )
  })

  crossed <- vapply(classifications, `[[`, logical(1L), "crossed")
  uncertain <- vapply(classifications, `[[`, logical(1L), "uncertain")
  current <- stage == "current"
  maximum <- stage == "maximum"
  unique_states <- sum(unique_state)
  repeated_states <- length(keys) - unique_states

  list(
    current_successes = sum(crossed[current]),
    maximum_successes = if (check_futility) {
      sum(crossed[maximum])
    } else {
      0L
    },
    inner_mc_uncertain_now = sum(uncertain[current]),
    inner_mc_uncertain_max = if (check_futility) {
      sum(uncertain[maximum])
    } else {
      0L
    },
    reuse = data.frame(
      method = method,
      bin_method = if (method == "bayes-bin") bin_method else NA_character_,
      cache_enabled = deterministic,
      analysis_requests = length(keys),
      unique_count_states = unique_states,
      repeated_count_states = repeated_states,
      repeated_state_rate = repeated_states / length(keys),
      cache_hits = if (deterministic) repeated_states else 0L,
      cache_hit_rate = if (deterministic) {
        repeated_states / length(keys)
      } else {
        0
      },
      stringsAsFactors = FALSE
    )
  )
}

#' Identify equivalent completed binary analyses
#'
#' @description Creates one identifier per predictive replicate from the arm-
#'   specific sufficient statistics and all fixed analysis settings. Replicates
#'   with the same identifier have the same deterministic completed-data
#'   result. The success threshold is applied afterwards and is not part of
#'   this result.
#'
#' @param states A data frame containing event counts and sample sizes for the
#'   control and treatment arms.
#' @inheritParams analyse_predictive_binary_counts
#'
#' @return A character vector with one analysis identifier per row of `states`.
#'
#' @keywords internal
#' @noRd
binary_count_analysis_keys <- function(
  states,
  method,
  single_arm,
  alternative,
  h0,
  prior_bin,
  bin_method,
  N_mcmc
) {
  numeric_key <- function(x) {
    paste(sprintf("%.17g", x), collapse = ",")
  }
  settings <- paste(
    method,
    single_arm,
    alternative,
    numeric_key(h0),
    numeric_key(prior_bin),
    bin_method,
    as.integer(N_mcmc),
    sep = "\r"
  )
  paste(
    settings,
    states$events_control,
    states$n_control,
    states$events_treatment,
    states$n_treatment,
    sep = "\r"
  )
}

#' Analyze one completed binary endpoint from sufficient statistics
#'
#' @param state A one-row data frame containing the event counts and sample
#'   sizes for the control and treatment arms.
#' @inheritParams analyse_predictive_binary_counts
#'
#' @return A list containing the completed-data success score and estimated
#'   treatment effect.
#'
#' @keywords internal
#' @noRd
analyse_binary_count_state_kernel <- function(
  state,
  single_arm,
  N_mcmc,
  method,
  alternative,
  h0,
  prior_bin,
  bin_method
) {
  if (method == "bayes-bin") {
    return(bayes_binomial_count_kernel(
      events_control = state$events_control,
      n_control = state$n_control,
      events_treatment = state$events_treatment,
      n_treatment = state$n_treatment,
      single_arm = single_arm,
      alternative = alternative,
      h0 = h0,
      prior_bin = prior_bin,
      bin_method = bin_method,
      N_mcmc = N_mcmc
    ))
  }

  fit <- risk_difference_wald_count_kernel(
    events_control = state$events_control,
    n_control = state$n_control,
    events_treatment = state$events_treatment,
    n_treatment = state$n_treatment,
    alternative = alternative,
    h0 = h0
  )
  list(success = fit$success, effect = fit$estimate)
}
