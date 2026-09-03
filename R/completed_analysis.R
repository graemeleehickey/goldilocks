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
    # The patient-level function computes the sufficient statistics once,
    # then uses the same prepared-data calculation as predictive imputation.
    # Keeping a single implementation prevents the two analysis routes from
    # drifting apart statistically.
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
    analysis <- analyse_logrank(
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
    analysis <- analyse_cox(
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

#' Analyze completed survival outcomes
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
analyse_completed_survival <- function(
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
    return(analyse_logrank(
      time = outcome$time,
      event = outcome$event,
      treatment = outcome$treatment,
      alternative = alternative
    ))
  }
  if (method == "cox") {
    return(analyse_cox(
      time = outcome$time,
      event = outcome$event,
      treatment = outcome$treatment,
      alternative = alternative,
      h0 = h0
    ))
  }
  if (method == "bayes-surv") {
    data_summ <- summarise_survival_outcomes(
      time = outcome$time,
      event = outcome$event,
      treatment = outcome$treatment,
      cutpoints = cutpoints,
      single_arm = single_arm
    )
    return(analyse_prepared_bayes_surv(
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
analyse_logrank <- function(time, event, treatment, alternative) {
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
#' @inheritParams analyse_logrank
#' @param h0 A finite numeric value giving the null log hazard ratio.
#'
#' @return A list containing `1 - P` for the Cox Wald test and the estimated log
#'   hazard ratio.
#'
#' @keywords internal
#' @noRd
analyse_cox <- function(time, event, treatment, alternative, h0) {
  fit_cox <- cox_wald_outcomes_checked(time, event, treatment)
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
