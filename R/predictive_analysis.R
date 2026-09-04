#' @title Evaluate success for one predictively completed trial
#'
#' @description Constructs one completed trial from the predictive imputations
#'   and evaluates the prespecified success criterion at the current and, when
#'   requested, maximum sample size.
#'
#' @inheritParams survival_adapt
#' @inheritParams sim_comp_data
#' @param prepared_outcomes A list of fixed outcomes and imputation
#'   positions returned by `prepare_predictive_outcomes()`.
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
#' @return A list containing completed-data analysis results and their success
#'   classifications at the current and, when requested, maximum sample sizes.
#' @noRd
analyse_predictive_survival <- function(
  prepared_outcomes,
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
  empty_interval,
  check_futility,
  prob_ha,
  mc_conf_level
) {
  ##############################################################################
  ### Test for success at current sample size (-> stop for success)
  ##############################################################################

  outcome_now <- complete_predictive_outcomes(
    prepared = prepared_outcomes,
    imputations = imputations,
    draw = draw
  )
  success_now <- analyse_completed_survival(
    outcome = outcome_now,
    cutpoints = cutpoints,
    end_of_study = end_of_study,
    interval_widths = interval_widths,
    prior_surv = prior_surv,
    N_mcmc = N_mcmc,
    single_arm = single_arm,
    method = method,
    alternative = alternative,
    h0 = h0,
    empty_interval = empty_interval
  )

  ##############################################################################
  ### Test for success at maximum sample size (-> stop for futility)
  ##############################################################################

  if (check_futility) {
    outcome_max <- complete_predictive_outcomes(
      prepared = prepared_outcomes,
      imputations = imputations,
      draw = draw,
      maximum = TRUE
    )
    success_max <- analyse_completed_survival(
      outcome = outcome_max,
      cutpoints = cutpoints,
      end_of_study = end_of_study,
      interval_widths = interval_widths,
      prior_surv = prior_surv,
      N_mcmc = N_mcmc,
      single_arm = single_arm,
      method = method,
      alternative = alternative,
      h0 = h0,
      empty_interval = empty_interval
    )
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
#' @param completed_counts A list returned by
#'   `predictive_binary_counts()`, containing arm-specific event counts
#'   and sample sizes for the current and maximum-sample completions.
#' @param single_arm A single logical value indicating a one-arm design.
#' @param N_mcmc A positive integer giving the number of beta-posterior draws
#'   for `method = "bayes-bin"` and `bin_method = "mc"`.
#' @param method A character value specifying `"riskdiff-wald"`,
#'   `"riskdiff-fm"`, or `"bayes-bin"`.
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
  completed_counts,
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
  n_draws <- nrow(completed_counts$current)
  if (check_futility) {
    combined <- rbind(completed_counts$current, completed_counts$maximum)
    draw_order <- as.vector(rbind(
      seq_len(n_draws),
      n_draws + seq_len(n_draws)
    ))
    combined <- combined[draw_order, , drop = FALSE]
    stage <- rep(c("current", "maximum"), n_draws)
  } else {
    combined <- completed_counts$current
    stage <- rep.int("current", n_draws)
  }
  rownames(combined) <- NULL

  analysis_groups <- group_binary_analyses(
    counts = combined,
    method = method,
    single_arm = single_arm,
    alternative = alternative,
    h0 = h0,
    prior_bin = prior_bin,
    bin_method = bin_method,
    N_mcmc = N_mcmc
  )
  deterministic <- method != "bayes-bin" || bin_method != "mc"
  first_in_group <- !duplicated(analysis_groups)
  analysis_rows <- if (deterministic) {
    which(first_in_group)
  } else {
    seq_len(nrow(combined))
  }

  analyses_unique <- lapply(analysis_rows, function(row) {
    analyse_binary_counts(
      counts = combined[row, , drop = FALSE],
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
    analyses_unique[match(
      analysis_groups,
      analysis_groups[first_in_group]
    )]
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
  unique_summaries <- sum(first_in_group)
  repeated_summaries <- length(analysis_groups) - unique_summaries

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
      reuse_enabled = deterministic,
      analysis_requests = length(analysis_groups),
      unique_count_summaries = unique_summaries,
      repeated_count_summaries = repeated_summaries,
      repeated_summary_rate = repeated_summaries / length(analysis_groups),
      reused_analyses = if (deterministic) repeated_summaries else 0L,
      reused_analysis_rate = if (deterministic) {
        repeated_summaries / length(analysis_groups)
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
#' @param counts A data frame containing event counts and sample sizes for the
#'   control and treatment arms.
#' @inheritParams analyse_predictive_binary_counts
#'
#' @return A character vector with one analysis identifier per row of `counts`.
#'
#' @keywords internal
#' @noRd
group_binary_analyses <- function(
  counts,
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
    counts$events_control,
    counts$n_control,
    counts$events_treatment,
    counts$n_treatment,
    sep = "\r"
  )
}

#' Analyze one completed binary endpoint from sufficient statistics
#'
#' @param counts A one-row data frame containing the event counts and sample
#'   sizes for the control and treatment arms.
#' @inheritParams analyse_predictive_binary_counts
#'
#' @return A list containing the completed-data success score and estimated
#'   treatment effect.
#'
#' @keywords internal
#' @noRd
analyse_binary_counts <- function(
  counts,
  single_arm,
  N_mcmc,
  method,
  alternative,
  h0,
  prior_bin,
  bin_method
) {
  if (method == "bayes-bin") {
    return(bayes_binomial_from_counts(
      events_control = counts$events_control,
      n_control = counts$n_control,
      events_treatment = counts$events_treatment,
      n_treatment = counts$n_treatment,
      single_arm = single_arm,
      alternative = alternative,
      h0 = h0,
      prior_bin = prior_bin,
      bin_method = bin_method,
      N_mcmc = N_mcmc
    ))
  }

  risk_difference_analysis <- if (method == "riskdiff-fm") {
    risk_difference_fm_from_counts
  } else {
    risk_difference_wald_from_counts
  }
  fit <- risk_difference_analysis(
    events_control = counts$events_control,
    n_control = counts$n_control,
    events_treatment = counts$events_treatment,
    n_treatment = counts$n_treatment,
    alternative = alternative,
    h0 = h0
  )
  list(success = fit$success, effect = fit$estimate)
}
