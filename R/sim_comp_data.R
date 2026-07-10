#' @title Simulate a complete clinical trial with event data drawn from a
#'   piecewise exponential distribution
#'
#' @param hazard_treatment vector. Finite non-negative constant hazard rates
#'   under the treatment arm.
#' @param hazard_control vector. Finite non-negative constant hazard rates
#'   under the control arm.
#' @param cutpoints finite, strictly increasing times at which the baseline
#'   hazard changes. The first value must be 0. Default is `cutpoints = 0`,
#'   which corresponds to a simple (non-piecewise) exponential model.
#' @param N_total integer. Maximum sample size allowable
#' @param lambda vector. Enrollment rates across simulated enrollment times. See
#'   [enrollment()] for more details.
#' @param lambda_time vector. Enrollment time(s) at which the enrollment rates
#'   change. Must be same length as lambda. See [enrollment()] for
#'   more details.
#' @param end_of_study finite study endpoint, strictly greater than the last
#'   cutpoint.
#' @param block scalar. Block size for generating the randomization schedule.
#' @param rand_ratio vector. Randomization allocation for the ratio of control
#'   to treatment. Integer values mapping the size of the block. See
#'   [randomization()] for more details.
#' @param prop_loss scalar. Overall proportion of subjects lost to follow-up.
#'   Subjects are selected at random for LTFU regardless of treatment assignment or
#'   event status. Each LTFU subject's observed time is drawn from a
#'   `Uniform(0, t)` distribution, where `t` is their potential
#'   event or censoring time. Since the LTFU time is always less than
#'   `t`, the event has not yet occurred at dropout and the subject is
#'   right-censored. Defaults to zero.
#'
#' @return A data frame with 1 row per subject and columns:
#'
#'   - `time`: Time of event or censoring time.
#'   - `treatment`: Treatment assignment, coded `1L` for the treatment arm and
#'     `0L` for the control arm. Single-arm designs have `treatment = 1L` for
#'     every subject.
#'   - `event`: Indicator of whether event occurred (`1L` if occurred and `0L`
#'     if right-censored).
#'   - `enrollment`: Time of patient enrollment relative to the time the trial
#'     enrolled the first patient. The package treats enrollment and
#'     randomization as occurring at the same time.
#'   - `id`: Identification number for each patient.
#'   - `loss_to_fu`: Indicator of whether the patient was lost to follow-up
#'     during observation.
#'
#' @importFrom stats runif sd
#' @export
sim_comp_data <- function(
  hazard_treatment,
  hazard_control = NULL,
  cutpoints = 0,
  N_total,
  lambda = 0.3,
  lambda_time = 0,
  end_of_study,
  block = 2,
  rand_ratio = c(1, 1),
  prop_loss = 0
) {
  ##############################################################################
  ### Run checks on arguments
  ##############################################################################

  validate_positive_integer_scalar(N_total, "N_total")
  validate_single_probability(prop_loss, "prop_loss")
  validate_cutpoints(cutpoints)
  validate_endpoint_time(end_of_study, cutpoints, "end_of_study")

  # Assign: indicator of whether single-arm study
  single_arm <- is.null(hazard_control)
  validate_piecewise_hazard(hazard_treatment, cutpoints, "hazard_treatment")
  if (!single_arm) {
    validate_piecewise_hazard(hazard_control, cutpoints, "hazard_control")
  }

  ##############################################################################
  ### Simulate enrollment/randomization + treatment assignment
  ##############################################################################

  # Simulate enrollment times
  enrollment <- enrollment(
    lambda = lambda,
    N_total = N_total,
    lambda_time = lambda_time
  )
  enrollment <- enrollment + runif(length(enrollment))
  enrollment <- sort(enrollment)

  # Simulate treatment assignment. The data convention is treatment = 1 for the
  # treatment arm and treatment = 0 for the control arm.
  if (!single_arm) {
    treatment <- randomization(
      N_total = N_total,
      block = block,
      allocation = rand_ratio
    )
  } else {
    treatment <- rep(1, N_total)
  }

  ##############################################################################
  ### Simulate event times
  ##############################################################################

  time <- rep(NA, length = N_total)
  event <- rep(NA, length = N_total)

  # Simulate TTE outcome
  # - Note: time = time *from* enrollment/randomization. In this package these
  #   are treated as the same time point.
  if (!single_arm) {
    sim_control <- pwe_sim(
      hazard = hazard_control,
      n = sum(treatment == 0),
      maxtime = end_of_study,
      cutpoints = cutpoints
    )
    time[treatment == 0] <- sim_control$time
    event[treatment == 0] <- sim_control$event
  }

  sim_treatment <- pwe_sim(
    hazard = hazard_treatment,
    n = sum(treatment == 1),
    maxtime = end_of_study,
    cutpoints = cutpoints
  )
  time[treatment == 1] <- sim_treatment$time
  event[treatment == 1] <- sim_treatment$event

  # Simulate loss to follow-up
  loss_to_fu <- rep(FALSE, N_total)
  if (prop_loss > 0) {
    n_loss_to_fu <- ceiling(prop_loss * N_total)
    loss_to_fu[sample(1:N_total, n_loss_to_fu)] <- TRUE
  }

  # Creating a new data.frame for all the variables
  data_total <- data.frame(
    time = time,
    treatment = treatment,
    event = event,
    enrollment = enrollment,
    id = 1:N_total,
    loss_to_fu = loss_to_fu
  )

  # Subjects lost are uniformly distributed
  if (prop_loss > 0) {
    data_total$time[data_total$loss_to_fu] <- runif(
      n_loss_to_fu,
      0,
      data_total$time[data_total$loss_to_fu]
    )
    data_total$event[data_total$loss_to_fu] <- rep(0, n_loss_to_fu)
  }

  return(data_total)
}
