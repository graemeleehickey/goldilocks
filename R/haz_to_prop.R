#' @title Posterior event estimates at endpoint time from posterior hazard
#'   distributions
#'
#' @description The posterior distribution of the probability of an event
#'   determined for a fixed time can be calculated in closed-form using the
#'   piecewise exponential distribution and a sample from the posterior
#'   distribution of the piecewise exponential hazard constant hazard rates.
#'
#' @inheritParams survival_adapt
#' @param post list. A list of posterior probabilities of the piecewise hazard
#'   rates (\eqn{\lambda_j}, for \code{j=1, \dots, J}) estimated by
#'   \code{\link{posterior}}. The list has two data frame elements:
#'   \code{post_treatment} and \code{post_control}. Each data frame is of length
#'   \code{N_mcmc} (see \code{\link{survival_adapt}}) and has \eqn{J} columns --
#'   one column per each constant hazard piece.
#' @param single_arm logical. If \code{TRUE}, trial is single arm. Else, if
#'   \code{FALSE}, it is a randomized two-arm trial.
#'
#' @return A data frame with 3 columns: the posterior probabilities of the event
#'   for the \code{treatment} arm, the posterior probabilities of the event for
#'   the \code{control} arm, and the posterior distribution of the treatment
#'   \code{effect}.
#'
#' @export
haz_to_prop <- function(post, cutpoint, end_of_study, single_arm) {

  control <- all(!is.na(post$post_control)) # Is there a control arm?

  if (length(cutpoint) == 1) {
    # Standard exponential for zero cutpoint
    p_treatment <- pexp(q = end_of_study,
                        rate = post$post_treatment)
    if (!single_arm) {
      p_control <- pexp(q = end_of_study,
                        rate = post$post_control)
    } else {
      p_control <- NA
    }
  } else {
    # PWE for >=1 break
    p_treatment <- bayesDP::ppexp(
      q = end_of_study,
      x = post$post_treatment,
      cuts = cutpoint)
    if (!single_arm) {
      p_control <- bayesDP::ppexp(
        q = end_of_study,
        x = post$post_control,
        cuts = cutpoint)
    } else {
      p_control <- NA
    }
  }

  if (!single_arm) {
    effect <- p_treatment - p_control
  } else {
    effect <- p_treatment
  }

  return(data.frame(p_treatment, p_control, effect))

}
