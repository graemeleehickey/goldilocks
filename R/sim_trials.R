#' @title Simulate one or more clinical trials subject to known design
#'   parameters and treatment effect
#'
#' @description Simulate multiple clinical trials with fixed input parameters,
#'   and tidly extract the relevant data to generate operating characteristics.
#'
#' @inheritParams survival_adapt
#' @param N_trials integer. Number of trials to simulate.
#'
#' @details This is basically a wrapper function for
#'   \code{\link{survival_adapt}}, whereby we repeatedly run the function for a
#'   independent number of trials (all with the same input design parameters and
#'   treatment effect).
#'
#' @return Data frame with 1 row per simulated trial and columns for key summary
#'   statistics.
#'
#' @export
#'
#' @examples
#' hc <- haz_est(c(0.20, 0.30), 12, 36)
#' ht <- haz_est(c(0.05, 0.15), 12, 36)
#'
#' out <- sim_trials(
#'   hazard_treatment = ht,
#'   hazard_control = hc,
#'   cutpoint = c(0, 12),
#'   N_total = 600,
#'   lambda = 20,
#'   lambda_time = NULL,
#'   interim_look = c(400, 500),
#'   end_of_study = 36,
#'   prior = c(.1, .1),
#'   block = 2,
#'   rand_ratio = c(1, 1),
#'   prop_loss_to_followup = 0.30,
#'   alternative = "less",
#'   h0 = 0,
#'   futility_prob = 0.05,
#'   expected_success_prob = 0.9,
#'   prob_ha = 0.975,
#'   N_impute = 5,
#'   N_mcmc = 5,
#'   N_trials = 2)

sim_trials <- function(
  hazard_treatment,
  hazard_control        = NULL,
  cutpoint              = 0,
  N_total,
  lambda                = 0.3,
  lambda_time           = NULL,
  interim_look          = NULL,
  end_of_study,
  prior                 = c(.1, .1),
  block                 = 2,
  rand_ratio            = c(1, 1),
  prop_loss_to_followup = 0.10,
  alternative           = "greater",
  h0                    = 0,
  futility_prob         = 0.05,
  expected_success_prob = 0.9,
  prob_ha               = 0.95,
  N_impute              = 10,
  N_mcmc                = 100,
  N_trials              = 10
) {

  # Setup empty vectors for output (length = number of simulated trials)
  N_enrolled            <- vector(length = N_trials)
  N_treatment           <- vector(length = N_trials)
  N_control             <- vector(length = N_trials)
  stop_futility         <- vector(length = N_trials)
  stop_expected_success <- vector(length = N_trials)
  est_interim           <- vector(length = N_trials)
  est_final             <- vector(length = N_trials)
  post_prob_ha          <- vector(length = N_trials)

  # Run adaptive design over simulated designs
  for (k in 1:N_trials) {
    sim <- survival_adapt(
      hazard_treatment      = hazard_treatment,
      hazard_control        = hazard_control,
      cutpoint              = cutpoint,
      N_total               = N_total,
      lambda                = lambda,
      lambda_time           = lambda_time,
      interim_look          = interim_look,
      end_of_study          = end_of_study,
      prior                 = prior,
      block                 = block,
      rand_ratio            = rand_ratio,
      prop_loss_to_followup = prop_loss_to_followup,
      alternative           = alternative,
      h0                    = h0,
      futility_prob         = futility_prob,
      expected_success_prob = expected_success_prob,
      prob_ha               = prob_ha,
      N_impute              = N_impute,
      N_mcmc                = N_mcmc
    )

    # Collect up the vectors
    N_enrolled[k]            <- sim$N_enrolled
    stop_futility[k]         <- sim$stop_futility
    stop_expected_success[k] <- sim$stop_expected_success
    est_interim[k]           <- sim$est_interim
    est_final[k]             <- sim$est_final
    post_prob_ha[k]          <- sim$post_prob_ha
  }

  out <- data.frame(
    N_enrolled,
    N_treatment,
    N_control,
    stop_futility,
    stop_expected_success,
    est_interim,
    est_final,
    post_prob_ha
  )

  return(out)

}
