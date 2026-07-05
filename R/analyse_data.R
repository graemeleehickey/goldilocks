#' @title Perform the final analysis test/method on the complete data
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
  bin_N = N_mcmc
) {
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
      single_arm = single_arm
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
    fit_cox <- cox_wald_test(data)
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

  treatment_stats <- beta_binomial_stats(
    event = data$event[data$treatment == 1],
    prior = bin_prior
  )

  if (single_arm) {
    if (bin_method == "mc") {
      effect <- rbeta(bin_N, treatment_stats$alpha, treatment_stats$beta)
      success <- posterior_tail_probability(effect, alternative, h0)
    } else {
      effect <- treatment_stats$mean
      success <- beta_binomial_single_arm_success(
        alpha = treatment_stats$alpha,
        beta = treatment_stats$beta,
        alternative = alternative,
        h0 = h0,
        method = bin_method
      )
    }
    return(list(success = success, effect = effect))
  }

  control_stats <- beta_binomial_stats(
    event = data$event[data$treatment == 0],
    prior = bin_prior
  )

  if (bin_method == "mc") {
    effect <- rbeta(bin_N, treatment_stats$alpha, treatment_stats$beta) -
      rbeta(bin_N, control_stats$alpha, control_stats$beta)
    success <- posterior_tail_probability(effect, alternative, h0)
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
    success <- beta_binomial_difference_success(
      treatment = treatment_stats,
      control = control_stats,
      alternative = alternative,
      h0 = h0
    )
  }

  list(success = success, effect = effect)
}

#' @noRd
validate_bayes_binomial_args <- function(bin_prior, bin_method, bin_N) {
  if (length(bin_prior) != 2 ||
    any(!is.finite(bin_prior)) ||
    any(bin_prior <= 0)) {
    stop("'bin_prior' must contain two positive finite values")
  }
  if (!bin_method %in% c("mc", "normal", "quadrature")) {
    stop("'bin_method' must be one of 'mc', 'normal', or 'quadrature'")
  }
  if (
    length(bin_N) != 1 ||
      !is.numeric(bin_N) ||
      is.na(bin_N) ||
      !is.finite(bin_N) ||
      bin_N <= 0 ||
      bin_N != floor(bin_N)
  ) {
    stop("'bin_N' must be a single positive integer")
  }
}

#' @noRd
beta_binomial_stats <- function(event, prior) {
  if (length(event) == 0) {
    stop("Bayesian binomial analysis requires at least one subject per arm")
  }
  alpha <- prior[1] + sum(event)
  beta <- prior[2] + length(event) - sum(event)
  list(
    alpha = alpha,
    beta = beta,
    mean = alpha / (alpha + beta),
    variance = (alpha * beta) / ((alpha + beta)^2 * (alpha + beta + 1))
  )
}

#' @noRd
posterior_tail_probability <- function(effect, alternative, h0) {
  if (alternative == "greater") {
    mean(effect > h0)
  } else if (alternative == "less") {
    mean(effect < h0)
  }
}

#' @noRd
beta_binomial_single_arm_success <- function(
  alpha,
  beta,
  alternative,
  h0,
  method
) {
  if (method == "normal") {
    effect <- alpha / (alpha + beta)
    effect_se <- sqrt((alpha * beta) / ((alpha + beta)^2 * (alpha + beta + 1)))
    if (alternative == "greater") {
      return(1 - pnorm(h0, mean = effect, sd = effect_se))
    } else if (alternative == "less") {
      return(pnorm(h0, mean = effect, sd = effect_se))
    }
  }

  if (alternative == "greater") {
    1 - pbeta(h0, alpha, beta)
  } else if (alternative == "less") {
    pbeta(h0, alpha, beta)
  }
}

#' @noRd
beta_binomial_difference_success <- function(
  treatment,
  control,
  alternative,
  h0
) {
  integrand <- function(x) {
    treatment_density <- dbeta(x, treatment$alpha, treatment$beta)
    control_threshold <- x - h0
    if (alternative == "greater") {
      treatment_density * pbeta(control_threshold, control$alpha, control$beta)
    } else if (alternative == "less") {
      treatment_density *
        (1 - pbeta(control_threshold, control$alpha, control$beta))
    }
  }

  integrate(integrand, lower = 0, upper = 1)$value
}

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

#' @title Fast Cox proportional hazards Wald test for treatment effect
#'
#' @inheritParams analyse_data
#'
#' @return A list with the estimated log hazard ratio (`estimate`) and its
#'   standard error (`std_error`).
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
