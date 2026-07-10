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
#'     - if `method = "chisq"`, 1 minus the chi-square test *P*-value.
#'   - `effect`: Sample vector from the posterior distribution of the effect
#'     size for `method = "bayes-surv"` and for `method = "bayes-bin"` with
#'     `bin_method = "mc"`, the posterior mean effect for `method =
#'     "bayes-bin"` with `bin_method = "normal"` or `"quadrature"`, the
#'     estimated log hazard ratio for `method = "cox"`, the chi-square
#'     statistic for `method = "chisq"`, or `NA` for `method = "logrank"`.
#'
#' @importFrom stats dbeta integrate pbeta pchisq pnorm rbeta
#' @import Rcpp
#' @import survival
#' @useDynLib goldilocks, .registration = TRUE
#'
#' @noRd
analyse_data <- function(
  data,
  cutpoints,
  end_of_study,
  prior,
  N_mcmc,
  single_arm,
  method,
  alternative,
  h0,
  bin_prior = c(1, 1),
  bin_method = "mc",
  bin_N = N_mcmc,
  empty_interval = "propagate"
) {
  validate_h0(h0, method, single_arm)

  ####################################################
  ### Bayesian test
  ####################################################

  # CIF_trt(T) - CIF_con(T) for two-armed trial

  if (method == "bayes-surv") {
    # Posterior distribution of lambdas: imputed data
    post_lambda_imp <- posterior(
      data = data,
      cutpoints = cutpoints,
      prior = prior,
      N_mcmc = N_mcmc,
      single_arm = single_arm,
      empty_interval = empty_interval
    )

    # Posterior distribution of event proportions: imputed data
    post_imp <- haz_to_prop(
      post = post_lambda_imp,
      cutpoints = cutpoints,
      end_of_study = end_of_study,
      single_arm = single_arm
    )

    effect <- post_imp$effect
    if (alternative == "greater") {
      success <- mean(effect > h0)
    } else if (alternative == "less") {
      success <- mean(effect < h0)
    }
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
      bin_prior = bin_prior,
      bin_method = bin_method,
      bin_N = bin_N
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
  ### Chi-square test
  ####################################################

  if (method == "chisq") {
    assert_complete_binary_outcomes(data, end_of_study, "Chi-square")
    mat <- with(data, table(event, treatment))
    fit_cs <- chisq.test(mat, correct = FALSE)
    success <- 1 - fit_cs$p.val
    effect <- fit_cs$statistic
  }

  return(list(
    "success" = success,
    "effect" = effect
  ))
}

#' @title Bayesian binomial test for complete binary outcomes
#'
#' @description Updates Beta priors for one or two treatment arms and calculates
#'   the posterior probability of the specified binary-outcome hypothesis.
#'
#' @inheritParams analyse_data
#'
#' @return A list with the posterior probability of success (`success`) and
#'   posterior treatment effect (`effect`).
#'
#' @noRd
bayes_binomial_test <- function(
  data,
  single_arm,
  alternative,
  h0,
  bin_prior,
  bin_method,
  bin_N
) {
  validate_bayes_binomial_args(bin_prior, bin_method, bin_N)
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
    alpha <- bin_prior[1] + sum(event)
    beta <- bin_prior[2] + length(event) - sum(event)
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
        treatment_density * pbeta(control_threshold, control$alpha, control$beta)
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
      effect <- rbeta(bin_N, treatment_stats$alpha, treatment_stats$beta)
      success <- if (alternative == "greater") {
        mean(effect > h0)
      } else {
        mean(effect < h0)
      }
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
    return(list(success = success, effect = effect))
  }

  control_stats <- beta_binomial_stats(data$event[data$treatment == 0])

  if (bin_method == "mc") {
    effect <- rbeta(bin_N, treatment_stats$alpha, treatment_stats$beta) -
      rbeta(bin_N, control_stats$alpha, control_stats$beta)
    success <- if (alternative == "greater") {
      mean(effect > h0)
    } else {
      mean(effect < h0)
    }
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

  list(success = success, effect = effect)
}

#' @title Validate complete binary outcomes
#'
#' @description Ensures a binary endpoint analysis is supplied only event
#'   indicators and that censored subjects have complete follow-up.
#'
#' @noRd
assert_complete_binary_outcomes <- function(data, end_of_study, method_label) {
  if (any(!data$event %in% c(0, 1))) {
    stop(method_label, " analysis requires binary event outcomes")
  }
  if (any(data$event == 0 & data$time < end_of_study)) {
    stop(
      method_label,
      " analysis requires all censored subjects to be followed ",
      "to 'end_of_study' or imputed before analysis"
    )
  }
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
  fit_warning <- NULL
  fit <- withCallingHandlers(
    cox_wald_test(data),
    warning = function(w) {
      fit_warning <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }
  )

  if (!is.null(fit_warning)) {
    stop(
      "Cox analysis is non-estimable: the Cox model did not produce a ",
      "reliable Wald statistic. survival::coxph.fit reported: ",
      fit_warning
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
