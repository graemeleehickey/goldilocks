#' @title Run test for whether we should stop enrollment at the interim analysis
#'
#' @description Imputes one predictive dataset and evaluates expected success
#'   at the current and, when requested, maximum sample size.
#'
#' @inheritParams survival_adapt
#' @inheritParams sim_comp_data
#' @param check_futility Logical. Does the adaptive design include a test for
#'   assessment of futility?
#' @param prob_ha Final completed-data success threshold.
#' @param mc_conf_level Confidence level for diagnostic finite posterior Monte
#'   Carlo bounds.
#'
#' @return See analyse_data
#' @noRd
test_stop_success <- function(
  data,
  hazard,
  end_of_study,
  cutpoints,
  single_arm,
  prior_surv,
  N_mcmc,
  method,
  alternative,
  h0,
  prior_bin,
  bin_method,
  binary_imputation,
  empty_interval,
  check_futility,
  prob_ha,
  mc_conf_level
) {
  ##############################################################################
  ### Test for success at current sample size (-> stop for success)
  ##############################################################################

  # Single imputed data set
  data_success_impute <- impute_data(
    data_in = data,
    hazard = hazard,
    end_of_study = end_of_study,
    cutpoints = cutpoints,
    type = "success",
    single_arm = single_arm,
    binary_imputation = binary_imputation
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
      empty_interval = empty_interval
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
    # Take the already imputed data for expected success and append on
    # imputed event times for subjects not yet enrolled

    # Single imputed data set
    data_futility_impute <- impute_data(
      data_in = data_success_impute,
      hazard = hazard,
      end_of_study = end_of_study,
      cutpoints = cutpoints,
      type = "futility",
      single_arm = single_arm,
      binary_imputation = binary_imputation
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
        empty_interval = empty_interval
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
