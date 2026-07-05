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
#'     - if `method = "bayes"`, the posterior probability that the treatment
#'       effect is greater than `h0` when `alternative = "greater"`, or less
#'       than `h0` when `alternative = "less"`;
#'     - if `method = "logrank"`, 1 minus the log-rank test *P*-value, using a
#'       two-sided *P*-value when `alternative = "two.sided"` and a one-sided
#'       *P*-value otherwise;
#'     - if `method = "cox"`, 1 minus the Cox Wald test *P*-value for the
#'       estimated log hazard ratio compared with `h0`, using a two-sided
#'       *P*-value when `alternative = "two.sided"` and a one-sided *P*-value
#'       otherwise;
#'     - if `method = "chisq"`, 1 minus the chi-square test *P*-value.
#'   - `effect`: Sample vector from the posterior distribution of the effect
#'     size for `method = "bayes"`, the estimated log hazard ratio for
#'     `method = "cox"`, the chi-square statistic for `method = "chisq"`, or
#'     `NA` for `method = "logrank"`.
#'
#' @importFrom stats pchisq pnorm
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
  h0
) {
  ####################################################
  ### Bayesian test
  ####################################################

  # CIF_trt(T) - CIF_con(T) for two-armed trial

  if (method == "bayes") {
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
    if (any(data$event == 0 & data$time < end_of_study)) {
      stop(
        "Chi-square analysis requires all censored subjects to be followed ",
        "to 'end_of_study' or imputed before analysis"
      )
    }
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
