#' @title Posterior event estimates at endpoint time from posterior hazard
#'   distributions
#'
#' @description The posterior distribution of the probability of an event
#'   determined for a fixed time can be calculated in closed-form using the
#'   piecewise exponential distribution and a sample from the posterior
#'   distribution of the piecewise exponential hazard constant hazard rates.
#'
#' @inheritParams survival_adapt
#' @inheritParams sim_comp_data
#' @param post array An array of posterior probabilities of the piecewise hazard
#'   rates (\eqn{\lambda_j}, for \code{j=1, \dots, J}) estimated by
#'   \code{\link{posterior}}. The array has dimension 3. The first dimension is
#'   of length \code{N_mcmc}, the second dimension is of length \eqn{J} (one
#'   column for each hazard piece), and the third dimension is of length 2, with
#'   the first slice including posterior samples from \code{post_treatment}, and
#'   the second slice including posterior samples from \code{post_control}.
#' @param single_arm logical. If \code{TRUE}, trial is single arm. Else, if
#'   \code{FALSE}, it is a randomized two-arm trial.
#'
#' @return A data frame with 3 columns of posterior samples:
#'
#'   \describe{
#'     \item{\code{p_treatment}}{
#'       The posterior probabilities of the event for the treatment arm.}
#'     \item{\code{p_control}}{
#'       The posterior probabilities of the event for the control arm.}
#'     \item{\code{effect}}{
#'       The posterior distribution of the treatment effect.}
#'   }
#'
#' @noRd
haz_to_prop <- function(post, cutpoints, end_of_study, single_arm) {

  if (length(cutpoints) == 1) {
    # Standard exponential for when no internal cutpoints
    p_treatment <- pexp(q    = end_of_study,
                        rate = post[, , 1])
    if (!single_arm) {
      p_control <- pexp(q    = end_of_study,
                        rate = post[, , 2])
    } else {
      p_control <- NA
    }
  } else {
    # PWE for >=1 break
    # p_treatment <- bayesDP::ppexp(
    #   q = end_of_study,
    #   x = post[, , 1],
    #   cuts = cutpoints)
    p_treatment <- ppwe(hazard = post[, , 1],
                        end_of_study = end_of_study,
                        cutpoints = cutpoints)
    if (!single_arm) {
      # p_control <- bayesDP::ppexp(
      #   q = end_of_study,
      #   x = post[, , 2],
      #   cuts = cutpoints)
      p_control <- ppwe(hazard = post[, , 2],
                        end_of_study = end_of_study,
                        cutpoints = cutpoints)
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
