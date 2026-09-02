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
