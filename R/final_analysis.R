#' @title Conduct the prespecified final analysis
#'
#' @description Applies the selected final analysis after accrual has stopped and
#'   the available follow-up is complete. Depending on `imputed_final`, subjects
#'   lost to follow-up are either handled as observed censored outcomes or have
#'   their outcomes multiply imputed.
#'
#' @inheritParams survival_adapt
#' @inheritParams sim_comp_data
#' @param data_in A data frame with one row per enrolled subject and columns for
#'   treatment assignment (`treatment`, coded `1`
#'   for treatment and `0` for control; single-arm designs use all 1s),
#'   event time (`time`), event indicator (`event`), and indicator of
#'   whether the subject requires imputation for expected success
#'   (`subject_impute_success`).
#'
#' @details Interim predictive calculations complete all outcomes that are not
#'   yet observed. At the final analysis, `imputed_final = TRUE` likewise
#'   imputes outcomes for subjects lost to follow-up before `end_of_study`.
#'   With `imputed_final = FALSE`, time-to-event analyses retain their observed
#'   right-censoring, whereas binary analyses (`method = "riskdiff-wald"`,
#'   `"riskdiff-fm"`, or `"bayes-bin"`) exclude subjects without complete
#'   endpoint status. Design evaluations should prespecify this choice and
#'   assess sensitivity to it when appreciable loss to follow-up is expected.
#'
#' @return A length-two numeric vector containing the posterior probability (or
#'   `1 - P` for a frequentist analysis) for the alternative hypothesis,
#'   followed by the treatment-effect estimate. For an imputed
#'   Cox analysis these are the Rubin-pooled Wald-test result and pooled log
#'   hazard ratio. For an imputed risk-difference analysis these are the
#'   Rubin-pooled Wald-test result and pooled treatment-control event-risk
#'   difference. Bayesian imputed analyses average their summaries over
#'   imputations.
#' @noRd
analyse_final <- function(
  data_in,
  cutpoints,
  prior_surv_final,
  N_mcmc,
  single_arm,
  imputed_final,
  method,
  N_impute,
  alternative,
  h0,
  prior_bin,
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
  interval_widths <- if (imputed_final && method == "bayes-surv") {
    endpoint_interval_widths(cutpoints, end_of_study)
  } else {
    NULL
  }

  if (imputed_final) {
    if (method %in% c("cox", "riskdiff-wald", "riskdiff-fm") && N_impute < 2) {
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
      prior_surv = prior_surv_final,
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

      if (method == "bayes-surv") {
        # The draw used by impute_data() came from the observed final-data
        # posterior. Preserve the documented two-stage procedure by forming a
        # new posterior from the completed imputation's sufficient statistics
        # and the final-stage prior. Passing statistics directly avoids rebuilding
        # a reduced patient-level analysis data frame for every imputation.
        data_summ <- posterior_sufficient_stats(
          data = data_success_impute,
          cutpoints = cutpoints,
          single_arm = single_arm
        )
        success <- analyse_prepared_bayes_surv(
          data_summ = data_summ,
          cutpoints = cutpoints,
          end_of_study = end_of_study,
          prior_surv = prior_surv_final,
          N_mcmc = N_mcmc,
          single_arm = single_arm,
          alternative = alternative,
          h0 = h0,
          empty_interval = empty_interval,
          interval_widths = interval_widths
        )
        post_paa[j] <- success$success
        effect_final[j] <- success$effect
        next
      }

      # The remaining methods consume subject-level outcomes or model fits.
      # Keep their established patient-level path unchanged.
      time <- NULL
      event <- NULL
      treatment <- NULL
      data <- subset(
        data_success_impute,
        select = c(time, event, treatment)
      )

      if (method == "cox") {
        fit_cox <- cox_wald_test_checked(data)
        effect_final[j] <- fit_cox$estimate
        variance_final[j] <- fit_cox$std_error^2
      } else if (method %in% c("riskdiff-wald", "riskdiff-fm")) {
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
          prior_surv = prior_surv_final,
          N_mcmc = N_mcmc,
          single_arm = single_arm,
          method = method,
          alternative = alternative,
          h0 = h0,
          prior_bin = prior_bin,
          bin_method = bin_method,
          empty_interval = empty_interval
        )

        post_paa[j] <- success$success
        if (method == "bayes-bin") {
          effect_final[j] <- success$effect
        }
      }
    }

    if (method %in% c("cox", "riskdiff-wald", "riskdiff-fm")) {
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
      method %in%
        c("riskdiff-wald", "riskdiff-fm", "bayes-bin") &&
        "loss_to_fu" %in% names(data_in)
    ) {
      data_in <- data_in[!data_in$loss_to_fu, ]
    }
    success <- analyse_data(
      data = data_in,
      cutpoints = cutpoints,
      end_of_study = end_of_study,
      prior_surv = prior_surv_final,
      N_mcmc = N_mcmc,
      single_arm = single_arm,
      method = method,
      alternative = alternative,
      h0 = h0,
      prior_bin = prior_bin,
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
#' @param estimates A numeric vector of treatment-effect estimates, one per
#'   imputed data set. At least two values are required.
#' @param variances A numeric vector of corresponding finite, non-negative
#'   within-imputation variances.
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
  if (!is.finite(total_variance) || total_variance < 0) {
    stop("Rubin pooling requires a finite non-negative total variance")
  }
  if (total_variance == 0) {
    difference_from_null <- estimate - h0
    statistic <- if (difference_from_null == 0) {
      0
    } else {
      sign(difference_from_null) * Inf
    }
    return(list(
      success = normal_test_success(statistic, alternative),
      estimate = estimate,
      std_error = 0,
      degrees_freedom = Inf
    ))
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
