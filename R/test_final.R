#' @title Run final test once enrollment is complete
#'
#' @description Once enrollment is complete, the final analysis must be
#'   conducted. The final analysis can take two forms depending on whether the
#'   data is "complete" or not (see Details).
#'
#' @inheritParams survival_adapt
#' @inheritParams sim_comp_data
#' @param data_in data frame. The time-to-event data that requires imputation,
#'   with columns for treatment assignment (`treatment`, coded `1`
#'   for treatment and `0` for control; single-arm designs use all 1s),
#'   event time (`time`), event indicator (`event`), and indicator of
#'   whether the subject requires imputation for expected success
#'   (`subject_impute_success`).
#'
#' @details If `prop_loss > 0`, then there is drop-out (loss to
#'   follow-up) for some subjects. During the Goldilocks algorithm, all subjects
#'   who are event-free and still eligible for follow-up for the primary
#'   endpoint are imputed. This means subjects who are lost to follow-up will be
#'   imputed. In the final analysis we need to decide whether we mimic that, or
#'   whether we will exclude or right-censor such subjects (depending on
#'   analysis `method` choice). This is controlled by `imputed_final`,
#'   which is a Boolean: `TRUE` implies the final analysis will be based on
#'   fully imputed data, whereas `FALSE` implies the final analysis will be
#'   based on observed (non-imputed) data only.
#'
#'   If `imputed_final = FALSE` then intuitively, we might expect to see a
#'   marginal increase in the proportion of studies that stop for expected
#'   success, but which then go on to fail. We have not verified this aspect,
#'   but it should be noted. When `method = "riskdiff"` or
#'   `method = "bayes-bin"` and `imputed_final = FALSE`, subjects lost to
#'   follow-up are excluded from the analysis because complete binary outcome
#'   analyses cannot handle right-censored observations.
#'
#' @return Vector with 1) posterior probability (or P-value equivalent) for
#'   the alternative hypothesis, and 2) the treatment effect. For an imputed
#'   Cox analysis these are the Rubin-pooled Wald-test result and pooled log
#'   hazard ratio. For an imputed risk-difference analysis these are the
#'   Rubin-pooled Wald-test result and pooled treatment-control event-risk
#'   difference. Bayesian imputed analyses average their summaries over
#'   imputations.
#' @noRd
test_final <- function(
  data_in,
  cutpoints,
  prior,
  N_mcmc,
  single_arm,
  imputed_final,
  method,
  N_impute,
  alternative,
  h0,
  bin_prior,
  bin_method,
  binary_imputation,
  empty_interval,
  end_of_study
) {
  validate_analysis_configuration(
    method,
    alternative,
    single_arm,
    imputed_final
  )

  if (imputed_final) {
    if (method %in% c("cox", "riskdiff") && N_impute < 2) {
      stop(
        "Frequentist final-analysis imputation requires at least two ",
        "imputations ",
        "to apply Rubin's rules"
      )
    }

    # Posterior distribution of lambdas: final data
    post_lambda_final <- posterior(
      data = data_in,
      cutpoints = cutpoints,
      prior = prior,
      N_mcmc = N_impute,
      single_arm = single_arm,
      empty_interval = empty_interval
    )
    # Effect estimate + posterior probability for each imputed dataset.
    # Frequentist model-based analyses additionally retain within-imputation
    # variances for Rubin pooling.
    effect_final <- rep(NA_real_, N_impute)
    post_paa <- rep(NA_real_, N_impute)
    variance_final <- rep(NA_real_, N_impute)
    # Impute multiple data sets
    for (j in 1:N_impute) {
      # Single imputed data set
      data_success_impute <- impute_data(
        data_in = data_in,
        hazard = post_lambda_final[j, , , drop = FALSE],
        end_of_study = end_of_study,
        cutpoints = cutpoints,
        type = "success",
        single_arm = single_arm,
        binary_imputation = binary_imputation
      )

      # Create enrolled subject data frame for analysis
      time <- NULL
      event <- NULL
      treatment <- NULL
      data <- subset(data_success_impute, select = c(time, event, treatment))

      if (method == "cox") {
        fit_cox <- cox_wald_test_checked(data)
        assert_cox_estimable(fit_cox)
        effect_final[j] <- fit_cox$estimate
        variance_final[j] <- fit_cox$std_error^2
      } else if (method == "riskdiff") {
        fit_riskdiff <- risk_difference_estimate_checked(
          data = data,
          end_of_study = end_of_study
        )
        effect_final[j] <- fit_riskdiff$estimate
        variance_final[j] <- fit_riskdiff$variance
      } else {
        # Apply primary analysis to imputed data
        success <- analyse_data(
          data = data,
          cutpoints = cutpoints,
          end_of_study = end_of_study,
          prior = prior,
          N_mcmc = N_mcmc,
          single_arm = single_arm,
          method = method,
          alternative = alternative,
          h0 = h0,
          bin_prior = bin_prior,
          bin_method = bin_method,
          empty_interval = empty_interval
        )

        post_paa[j] <- success$success
        if (method %in% c("bayes-surv", "bayes-bin")) {
          effect_final[j] <- success$effect
        }
      }
    }

    if (method %in% c("cox", "riskdiff")) {
      pooled <- pool_rubin_scalar(
        estimates = effect_final,
        variances = variance_final,
        alternative = alternative,
        h0 = h0
      )
      post_paa <- pooled$success
      est_final <- pooled$estimate
    } else {
      # Average Bayesian summaries over imputations
      post_paa <- mean(post_paa)
      est_final <- mean(effect_final)
    }
  } else {
    # Apply primary analysis to final data (without imputation)
    # Risk-difference and Bayesian binomial analyses cannot handle censored
    # (LTFU) subjects, so exclude them.
    if (
      method %in% c("riskdiff", "bayes-bin") && "loss_to_fu" %in% names(data_in)
    ) {
      data_in <- data_in[!data_in$loss_to_fu, ]
    }
    success <- analyse_data(
      data = data_in,
      cutpoints = cutpoints,
      end_of_study = end_of_study,
      prior = prior,
      N_mcmc = N_mcmc,
      single_arm = single_arm,
      method = method,
      alternative = alternative,
      h0 = h0,
      bin_prior = bin_prior,
      bin_method = bin_method,
      empty_interval = empty_interval
    )

    post_paa <- success$success
    est_final <- success$effect
  }

  return(c(post_paa, est_final))
}

#' @title Pool scalar treatment effects using Rubin's rules
#'
#' @description Combines scalar treatment-effect estimates and their
#'   within-imputation variances, then evaluates the pooled estimate against
#'   its null value using Rubin's large-sample degrees of freedom.
#'
#' @param estimates Numeric vector of per-imputation effect estimates.
#' @param variances Numeric vector of corresponding within-imputation variances.
#' @inheritParams survival_adapt
#'
#' @return A list containing the pooled success score, effect estimate, standard
#'   error, and degrees of freedom.
#'
#' @importFrom stats pt var
#' @noRd
pool_rubin_scalar <- function(estimates, variances, alternative, h0) {
  m <- length(estimates)
  if (m < 2 || length(variances) != m) {
    stop("Rubin pooling requires at least two paired estimates and variances")
  }
  if (
    anyNA(estimates) ||
      anyNA(variances) ||
      any(!is.finite(estimates)) ||
      any(!is.finite(variances)) ||
      any(variances < 0)
  ) {
    stop(
      "Rubin pooling requires finite estimates and non-negative variances"
    )
  }

  estimate <- mean(estimates)
  within_variance <- mean(variances)
  between_variance <- var(estimates)
  total_variance <- within_variance + (1 + 1 / m) * between_variance
  if (!is.finite(total_variance) || total_variance <= 0) {
    stop("Rubin pooling requires a finite positive total variance")
  }

  relative_increase <- if (within_variance == 0) {
    Inf
  } else {
    (1 + 1 / m) * between_variance / within_variance
  }
  degrees_freedom <- if (relative_increase == 0) {
    Inf
  } else {
    (m - 1) * (1 + 1 / relative_increase)^2
  }

  statistic <- (estimate - h0) / sqrt(total_variance)
  success <- switch(
    alternative,
    "less" = 1 - pt(statistic, df = degrees_freedom),
    "greater" = pt(statistic, df = degrees_freedom),
    "two.sided" = 1 - 2 * pt(-abs(statistic), df = degrees_freedom)
  )

  list(
    success = success,
    estimate = estimate,
    std_error = sqrt(total_variance),
    degrees_freedom = degrees_freedom
  )
}
