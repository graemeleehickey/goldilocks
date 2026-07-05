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
#' | Element | Description |
#' | --- | --- |
#' | `success` | Mean posterior probability of effect, or 1 minus the conventional *P*-value if `method = "logrank"`, `method = "cox"`, or `method = "chisq"`. |
#' | `effect` | Sample vector from the posterior distribution of the effect size. If `method = "logrank"`, this is `NULL`. |
#'
#' @importFrom stats pchisq pnorm
#' @import survival
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
    fit_cox <- coxph(Surv(time, event) ~ treatment, data = data)
    fit_res <- coef(summary(fit_cox))
    z <- (fit_res[1, 1] - h0) / fit_res[1, 3]
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
    effect <- fit_res[1, 1]
  }

  ####################################################
  ### Chi-square test
  ####################################################

  # Assumes all LTFU subjects have been imputed

  if (method == "chisq") {
    mat <- with(data, table(event, treatment))
    fit_cs <- chisq.test(mat, correct = FALSE)
    success <- 1 - fit_cs$p.val
    effect <- fit_cs$statistic
  }

  return(list(
    "success" = success,
    "effect"  = effect
  ))
}
