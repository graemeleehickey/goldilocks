#' @title Convert posterior hazard draws to endpoint event probabilities
#'
#' @description Calculates arm-specific event probabilities at a fixed follow-up
#'   time from draws of piecewise-exponential hazard rates, together with the
#'   treatment effect on the event-probability scale.
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
#'     treatment-minus-control event-probability difference for a two-arm
#'     design.
#'
#' @noRd
haz_to_prop <- function(post, cutpoints, end_of_study, single_arm) {
  validate_cutpoints(cutpoints)
  validate_endpoint_time(end_of_study, cutpoints, "end_of_study")
  hazards <- bayes_surv_cum_haz_draws(
    post, endpoint_interval_widths(cutpoints, end_of_study), single_arm
  )
  p_treatment <- cum_haz_to_prob(hazards$treatment)
  p_control <- if (single_arm) {
    NA_real_
  } else {
    cum_haz_to_prob(hazards$control)
  }
  # Calculate the difference before rounding arm risks for display. Both risks
  # can be reported as one while their difference remains nonzero.
  effect <- if (single_arm) {
    p_treatment
  } else {
    bayes_surv_risk_difference(hazards$treatment, hazards$control)$effect
  }
  data.frame(p_treatment, p_control, effect)
}
