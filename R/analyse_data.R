#' @title Perform the final analysis test/method on the complete data
#'
#' @description Dispatches an imputed or complete trial dataset to the selected
#'   Bayesian or frequentist analysis and returns its success score and effect.
#'
#' @inheritParams survival_adapt
#' @inheritParams sim_comp_data
#' @inheritParams haz_to_prop
#' @param data data frame. The (time-to-event) analysis data to be analyzed per
#'   the pre-specified analysis method. Generally this will be an imputed data
#'   set, and the analysis will be looped over multiple imputed datasets.
#'
#' @return A list with 2 elements:
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
    t0 <- data$time[data$treatment == 0]
    t1 <- data$time[data$treatment == 1]
    e0 <- data$event[data$treatment == 0]
    e1 <- data$event[data$treatment == 1]
    lr <- logrank_test(t0, t1, e0, e1)
    assert_logrank_estimable(lr)
    if (alternative == "two.sided") {
      success <- 1 - lr[3]
    } else {
      # Logrank z > 0 when control has excess events (treatment beneficial).
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
    effect <- NA
  }

  ####################################################
  ### Cox regression test
  ####################################################

  if (method == "cox") {
    fit_cox <- cox_wald_test_checked(data)
    assert_cox_estimable(fit_cox)
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
    effect <- fit_cox$estimate
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

#' @title Analyze Bayesian survival sufficient statistics
#'
#' @description Applies the Bayesian piecewise-exponential final analysis
#'   directly to completed-data sufficient statistics. Predictive-imputation
#'   callers use this function after deriving the completed exposure and event
#'   totals, avoiding a second patient-level posterior setup.
#'
#' @param data_summ Completed-data sufficient statistics created by
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

  # Keep the posterior arithmetic next to the analysis that consumes it; the
  # two helpers capture the only non-trivial repeated calculations.
  beta_binomial_stats <- function(event) {
    if (length(event) == 0) {
      stop("Bayesian binomial analysis requires at least one subject per arm")
    }
    alpha <- prior_bin[1] + sum(event)
    beta <- prior_bin[2] + length(event) - sum(event)
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

  treatment_stats <- beta_binomial_stats(data$event[data$treatment == 1])

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

  control_stats <- beta_binomial_stats(data$event[data$treatment == 0])

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

  treatment_risk <- mean(treatment_event)
  control_risk <- mean(control_event)
  estimate <- treatment_risk - control_risk
  variance <- treatment_risk *
    (1 - treatment_risk) /
    length(treatment_event) +
    control_risk * (1 - control_risk) / length(control_event)

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

#' @title Calculate the log-rank test very quickly
#'
#' @description Calls the compiled log-rank implementation for two groups of
#'   survival outcomes and event indicators.
#'
#' @param groupa vector of group a's survival times
#' @param groupb vector of group b's survival times
#' @param groupacensored vector of censored information of group a's survival
#'   times
#' @param groupbcensored vector of censored information of group b's survival
#'   times
#' @param onlyz (optional) calculate only z-statistic
#'
#' @return chi-squared statistic, z-statistic, p-value
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

#' @title Fast Cox proportional hazards Wald test for treatment effect
#'
#' @description Fits the package's low-overhead Cox model and converts fitter
#'   warnings into clear non-estimability errors.
#'
#' @inheritParams analyse_data
#'
#' @return A list with the estimated log hazard ratio (`estimate`) and its
#'   standard error (`std_error`).
#'
#' @noRd
cox_wald_test_checked <- function(data) {
  fit_state <- new.env(parent = emptyenv())
  fit_state$warning <- NULL
  fit <- withCallingHandlers(
    cox_wald_test(data),
    warning = function(w) {
      fit_state$warning <- conditionMessage(w)
      invokeRestart("muffleWarning")
    }
  )

  if (!is.null(fit_state$warning)) {
    stop(
      "Cox analysis is non-estimable: the Cox model did not produce a ",
      "reliable Wald statistic. survival::coxph.fit reported: ",
      fit_state$warning
    )
  }

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

#' @title Fit a low-overhead Cox Wald test
#'
#' @description Calls survival's internal Cox fitter directly and returns the
#'   treatment log hazard ratio with its standard error.
#'
#' @noRd
cox_wald_test <- function(data) {
  y <- Surv(data$time, data$event)
  x <- matrix(as.double(data$treatment), ncol = 1)

  # Use the lower-level fitter to avoid formula/model-frame and summary
  # overhead in repeated simulation analyses.
  coxph_fit <- getFromNamespace("coxph.fit", "survival")
  fit <- coxph_fit(
    x = x,
    y = y,
    strata = NULL,
    offset = NULL,
    init = NULL,
    control = coxph.control(),
    weights = NULL,
    method = "efron",
    rownames = NULL,
    resid = FALSE,
    nocenter = NULL
  )

  list(
    estimate = fit$coefficients[1],
    std_error = sqrt(fit$var[1, 1])
  )
}
