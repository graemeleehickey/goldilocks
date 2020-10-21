#' @title Perform the final analysis test/method on the complete data
#'
#' @inheritParams survival_adapt
#' @inheritParams haz_to_prop
#' @param data data frame. The time-to-event analysis data to be analyzed per
#'   the pre-specified analysis method. Generally this will be an imputed data
#'   set, and the analysis will be looped over multiple imputed datasets.
#'
#' @return A list with 2 elements:
#'
#' \describe{
#'     \item{\code{success}}{
#'        The mean posterior probability of effect (or 1 - the conventional
#'        P-value if \code{method = "logrank"}).}
#'     \item{\code{effect}}{
#'        A sample from the posterior distribution of the effect size. If
#'        \code{method = "logrank"}, then this is \code{NULL}.}
#' }
#'
#' @importFrom stats pchisq
#' @export
analyse_data <- function(data, cutpoint, prior, N_mcmc, single_arm, method) {

  # If we do a Bayesian test of CIF_mesh(T) - CIF_control(T) on each
  # imputed data set
  if (method == "bayes") {
    # Posterior distribution of lambdas: imputed data
    post_lambda_imp <- posterior(data       = data,
                                 cutpoint   = cutpoint,
                                 prior      = prior,
                                 N_mcmc     = N_mcmc,
                                 single_arm = single_arm)

    # Posterior distribution of event proportions: imputed data
    post_imp <- haz_to_prop(post         = post_lambda_imp,
                            cutpoint     = cutpoint,
                            end_of_study = end_of_study,
                            single_arm   = single_arm)

    # Apply statistical test to declare success (e.g. efficacy)
    effect <- post_imp$effect

    # Apply test with direction
    if (alternative == "two-sided") {
      success <- max(c(mean(effect > h0), mean(effect < h0)))
    } else if (alternative == "greater") {
      success <- mean(effect > h0)
    } else {
      success <- mean(effect < h0)
    }

  }

  # If we do a log-rank test on each imputed data set
  if (method == "logrank") {
    lrt <- survdiff(Surv(time, event) ~ treatment, data = data)
    p <- pchisq(lrt$chisq, length(lrt$n) - 1, lower.tail = FALSE)
    success <- 1 - p
    effect <- NULL
  }

  return(list(
    "success" = success,
    "effect"  = effect))

}
