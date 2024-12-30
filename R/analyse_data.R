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
#' \describe{
#'     \item{\code{success}}{
#'        The mean posterior probability of effect (or 1 - the conventional
#'        P-value if \code{method = "logrank"}, \code{method = "cox"}, or
#'        \code{method = "chisq"}).}
#'     \item{\code{effect}}{ A sample vector from the posterior distribution of
#'        the effect size. If \code{method = "logrank"}, then this is
#'        \code{NULL}.}
#' }
#'
#' @importFrom stats pchisq
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
  h0) {

  ####################################################
  ### Bayesian test
  ####################################################

  # CIF_trt(T) - CIF_con(T) for two-armed trial

  if (method == "bayes") {
    # Posterior distribution of lambdas: imputed data
    post_lambda_imp <- posterior(data       = data,
                                 cutpoints  = cutpoints,
                                 prior      = prior,
                                 N_mcmc     = N_mcmc,
                                 single_arm = single_arm)

    # Posterior distribution of event proportions: imputed data
    post_imp <- haz_to_prop(post         = post_lambda_imp,
                            cutpoints    = cutpoints,
                            end_of_study = end_of_study,
                            single_arm   = single_arm)

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
    p  <- logrank_test(t0, t1, e0, e1)[3]
    # lrt <- survdiff(Surv(time, event) ~ treatment, data = data)
    # p <- pchisq(lrt$chisq, 1, lower.tail = FALSE)
    success <- 1 - p
    effect <- NA
  }

  ####################################################
  ### Cox regression test
  ####################################################

  if (method == "cox") {
    fit_cox <- coxph(Surv(time, event) ~ treatment, data = data)
    fit_res <- coef(summary(fit_cox))
    success <- 1 - fit_res[5]
    effect <- fit_res[1]
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
    "effect"  = effect))

}
