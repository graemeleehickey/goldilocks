#' @title Convert posterior hazard draws to endpoint event probabilities
#'
#' @description Calculates arm-specific event probabilities at a fixed
#'   follow-up time from draws of piecewise-exponential hazard rates, together
#'   with the treatment effect on the event-probability scale.
#'
#' @inheritParams survival_adapt
#' @inheritParams sim_comp_data
#' @param post A three-dimensional numeric array of posterior hazard-rate draws
#'   from `posterior()`. The dimensions are posterior draws, piecewise
#'   intervals, and arms, with treatment in the first arm slice and control in
#'   the second.
#' @param single_arm A single logical value indicating whether the trial has one
#'   treatment arm (`TRUE`) or treatment and control arms (`FALSE`).
#'
#' @return A data frame with three columns of posterior draws:
#'
#'   - `p_treatment`: Posterior probabilities of the event for the treatment
#'     arm.
#'   - `p_control`: Posterior event probabilities for the control arm, or `NA`
#'     for a single-arm design.
#'   - `effect`: Treatment-arm event probability for a single-arm design, or the
#'     treatment-minus-control event-probability difference for a two-arm design.
#'
#' @noRd
haz_to_prop <- function(post, cutpoints, end_of_study, single_arm) {
  if (length(cutpoints) == 0) {
    # Standard exponential for when there are no interior cutpoints
    p_treatment <- cumulative_hazard_to_probability(
      end_of_study * post[,, 1]
    )
    if (!single_arm) {
      p_control <- cumulative_hazard_to_probability(
        end_of_study * post[,, 2]
      )
    } else {
      p_control <- NA
    }
  } else {
    # Piecewise exponential for one or more interior cutpoints
    # Preserve the posterior-draw dimension when N_mcmc = 1. Direct array
    # slicing drops it to a vector, whereas ppwe() requires a hazard matrix.
    treatment_hazard <- matrix(
      post[,, 1],
      nrow = dim(post)[1],
      ncol = dim(post)[2]
    )
    p_treatment <- ppwe(
      hazard = treatment_hazard,
      end_of_study = end_of_study,
      cutpoints = cutpoints
    )
    if (!single_arm) {
      control_hazard <- matrix(
        post[,, 2],
        nrow = dim(post)[1],
        ncol = dim(post)[2]
      )
      p_control <- ppwe(
        hazard = control_hazard,
        end_of_study = end_of_study,
        cutpoints = cutpoints
      )
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
