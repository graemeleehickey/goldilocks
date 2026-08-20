#' @title Simulate a complete clinical trial with event data drawn from a
#'   piecewise exponential distribution
#'
#' @param hazard_treatment vector. Finite non-negative constant hazard rates
#'   under the treatment arm.
#' @param hazard_control vector. Finite non-negative constant hazard rates
#'   under the control arm.
#' @param cutpoints finite, positive, strictly increasing interior times at
#'   which the baseline hazard changes. The number of hazards for each arm must
#'   be one greater than the number of cutpoints. Default is `NULL`, which
#'   corresponds to a simple (non-piecewise) exponential model. Realized event
#'   times are assigned to analysis intervals using the survival
#'   counting-process convention, open on the left and closed on the right; an
#'   event exactly at a cutpoint therefore belongs to the interval ending at
#'   that cutpoint.
#' @param N_total integer. Maximum sample size allowable
#' @param lambda finite positive enrollment rates per unit time. Supply one rate
#'   for each interval defined by `lambda_time`. See [enrollment()] for the
#'   precise continuous-time process and time-origin convention.
#' @param lambda_time `NULL`, or finite, positive, strictly increasing internal
#'   times at which the enrollment rate changes. The initial boundary at zero
#'   is implicit, so `length(lambda)` must equal `length(lambda_time) + 1`.
#' @param end_of_study finite study endpoint, strictly greater than the last
#'   cutpoint.
#' @param block scalar. Block size for generating the randomization schedule.
#' @param rand_ratio length-two positive integer randomization allocation. Name
#'   the values `control` and `treatment`; either supplied order is accepted and
#'   normalized internally. A legacy unnamed vector remains accepted in
#'   `c(control, treatment)` order. Unequal unnamed values produce a warning
#'   because names may be required in a future major release. See
#'   [randomization()] for more details.
#' @param prop_loss one or two probabilities. A scalar applies the same LTFU
#'   proportion to every arm. For a two-arm design, differential attrition can
#'   be specified with a length-two vector named `control` and `treatment`; the
#'   supplied order does not matter. Within each arm,
#'   `ceiling(prop_loss * arm size)` subjects are selected at random regardless
#'   of event status. Each selected subject's observed time is drawn from a
#'   `Uniform(0, t)` distribution, where `t` is their potential event or
#'   censoring time. Since the LTFU time is always less than `t`, the event has
#'   not yet occurred at dropout and the subject is right-censored. Single-arm
#'   designs require one probability. Defaults to zero.
#'
#' @details Enrollment is simulated directly in continuous time by
#'   [enrollment()]. The first patient is placed at time zero and all subsequent
#'   enrollment times are measured from first patient in. No uniform jitter is
#'   added in `sim_comp_data()`.
#'
#'   `lambda_time` and `cutpoints` both contain internal change times, but they
#'   describe different clocks. `lambda_time` describes changes in the trial's
#'   calendar-time enrollment rate measured from first patient in. `cutpoints`
#'   describes changes in an individual subject's event hazard measured from
#'   that subject's enrollment. They need not have the same values or length.
#'   All time quantities supplied to a simulation should nevertheless use one
#'   common unit, such as days or months.
#'
#'   PWEALL represents the continuous generating hazard with pieces closed on
#'   the left and open on the right. This differs from the package's
#'   open-left, closed-right convention for assigning realized times only at
#'   the cutpoints themselves, which have probability zero under the continuous
#'   model. The cumulative hazard, event-time distribution, and generated
#'   simulations are therefore unchanged.
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
  cutpoints = NULL,
  N_total,
  lambda = 0.3,
  lambda_time = NULL,
  end_of_study,
  block = 2,
  rand_ratio = c(control = 1, treatment = 1),
  prop_loss = 0
) {
  ##############################################################################
  ### Run checks on arguments
  ##############################################################################

  # Assign: indicator of whether single-arm study
  single_arm <- is.null(hazard_control)

  validate_positive_integer_scalar(N_total, "N_total")
  prop_loss <- normalize_prop_loss(prop_loss, single_arm)
  validate_cutpoints(cutpoints)
  validate_endpoint_time(end_of_study, cutpoints, "end_of_study")

  validate_piecewise_hazard(hazard_treatment, cutpoints, "hazard_treatment")
  if (!single_arm) {
    validate_piecewise_hazard(hazard_control, cutpoints, "hazard_control")
    rand_ratio <- validate_randomization_args(
      N_total,
      block,
      rand_ratio,
      allocation_name = "rand_ratio"
    )
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
  treatment_values <- c(control = 0L, treatment = 1L)
  for (arm in names(prop_loss)) {
    arm_index <- which(treatment == treatment_values[[arm]])
    n_arm_loss <- ceiling(prop_loss[[arm]] * length(arm_index))
    if (n_arm_loss > 0L) {
      loss_to_fu[sample(arm_index, n_arm_loss)] <- TRUE
    }
  }
  n_loss_to_fu <- sum(loss_to_fu)

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
  if (n_loss_to_fu > 0L) {
    data_total$time[data_total$loss_to_fu] <- runif(
      n_loss_to_fu,
      0,
      data_total$time[data_total$loss_to_fu]
    )
    data_total$event[data_total$loss_to_fu] <- rep(0, n_loss_to_fu)
  }

  return(data_total)
}
