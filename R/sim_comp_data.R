#' @title Simulate complete trial data under piecewise-exponential event rates
#'
#' @description Simulates enrollment, treatment allocation, event or censoring
#'   times, and loss to follow-up for a single-arm or randomized two-arm trial.
#'   Event times follow a piecewise-exponential distribution within each arm.
#'
#' @param hazard_treatment A required numeric vector of finite, non-negative
#'   event rates for the treatment arm. Supply one rate per interval defined by
#'   `generation_cutpoints`; a single value specifies a constant event rate.
#' @param hazard_control `NULL` (the default) for a single-arm trial, or a
#'   numeric vector of finite, non-negative event rates for the control arm in
#'   a two-arm trial. It must contain one rate per interval defined by
#'   `generation_cutpoints`.
#' @param generation_cutpoints `NULL` (the default), or a numeric vector of
#'   finite, positive, strictly increasing interior
#'   follow-up times at which the data-generating hazard changes. The number of
#'   hazards for each arm must be one greater than the number of generation
#'   cutpoints. `NULL` specifies a constant-hazard data-generating model.
#' @param N_total A required positive integer giving the maximum total sample
#'   size.
#' @param lambda A numeric vector of finite, positive enrollment rates per unit
#'   of calendar time. Supply one rate for each interval defined by
#'   `lambda_time`. The default is `0.3`. See [enrollment()] for the
#'   continuous-time enrollment model and time origin.
#' @param lambda_time `NULL` (the default), or a numeric vector of finite,
#'   positive, strictly increasing calendar times at which the enrollment rate
#'   changes. Time zero is implicit, and `length(lambda)` must equal
#'   `length(lambda_time) + 1`.
#' @param end_of_study A required finite, positive numeric value giving the
#'   planned follow-up time for each subject. It must be later than the final
#'   `generation_cutpoints` value and use the same time unit.
#' @param block A positive integer vector of permitted randomization block
#'   sizes. Every value must be a multiple of `sum(rand_ratio)`. The default is
#'   `2` and the argument is ignored for a single-arm trial.
#' @param rand_ratio A length-two positive integer vector giving the control to
#'   treatment randomization ratio. The default is
#'   `c(control = 1, treatment = 1)`. Name
#'   the values `control` and `treatment`; either supplied order is accepted and
#'   matched by name. A legacy unnamed vector remains accepted in
#'   `c(control, treatment)` order. Unequal unnamed values produce a warning
#'   because names may be required in a future major release. See
#'   [randomization()] for more details.
#' @param prop_loss A numeric vector containing one or two probabilities in
#'   `[0, 1]`. A single value applies the same loss-to-follow-up proportion to
#'   every arm. For a two-arm design, differential attrition can
#'   be specified with a length-two vector named `control` and `treatment`; the
#'   supplied order does not matter. Within each arm,
#'   `ceiling(prop_loss * arm size)` subjects are selected at random regardless
#'   of event status. Each selected subject's observed time is drawn from a
#'   `Uniform(0, t)` distribution, where `t` is their potential event or
#'   censoring time. Since the LTFU time is always less than `t`, the event has
#'   not yet occurred at dropout and the subject is right-censored. Single-arm
#'   designs require one probability. The default is `0`, denoting no loss to
#'   follow-up.
#'
#' @details Enrollment is simulated directly in continuous time by
#'   [enrollment()]. The first patient is placed at time zero and all subsequent
#'   enrollment times are measured from first patient in. No uniform jitter is
#'   added in `sim_comp_data()`.
#'
#'   `lambda_time` and `generation_cutpoints` both contain internal change
#'   times, but they describe different clocks. `lambda_time` describes changes
#'   in the trial's calendar-time enrollment rate measured from first patient
#'   in. `generation_cutpoints` describes changes in an individual subject's
#'   event hazard measured from that subject's enrollment. They need not have
#'   the same values or length. All time quantities supplied to a simulation
#'   should nevertheless use one common unit, such as days or months.
#'
#'   PWEALL represents the continuous generating hazard with pieces closed on
#'   the left and open on the right. This differs from the package's
#'   open-left, closed-right convention for assigning realized times only at
#'   the cutpoints themselves, which have probability zero under the continuous
#'   model. The cumulative hazard, event-time distribution, and generated
#'   simulations are therefore unchanged.
#'
#' @return A data frame with one row per subject and columns:
#'
#'   - `time`: Numeric event or censoring time.
#'   - `treatment`: Numeric treatment indicator, coded `1` for the treatment arm
#'     and `0` for the control arm. Single-arm designs have `treatment = 1` for
#'     every subject.
#'   - `event`: Numeric event indicator, coded `1` for an event and `0` for
#'     right-censoring.
#'   - `enrollment`: Numeric time of subject enrollment relative to first
#'     patient in. The package treats enrollment and
#'     randomization as occurring at the same time.
#'   - `id`: Integer subject identifier.
#'   - `loss_to_fu`: Logical indicator of loss to follow-up.
#'
#' @importFrom stats runif sd
#' @export
sim_comp_data <- function(
  hazard_treatment,
  hazard_control = NULL,
  generation_cutpoints = NULL,
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
  validate_cutpoints(generation_cutpoints, "generation_cutpoints")
  validate_endpoint_time(
    end_of_study,
    generation_cutpoints,
    "end_of_study",
    "generation_cutpoints"
  )

  validate_piecewise_hazard(
    hazard_treatment,
    generation_cutpoints,
    "hazard_treatment",
    "generation_cutpoints"
  )
  if (!single_arm) {
    validate_piecewise_hazard(
      hazard_control,
      generation_cutpoints,
      "hazard_control",
      "generation_cutpoints"
    )
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
      cutpoints = generation_cutpoints
    )
    time[treatment == 0] <- sim_control$time
    event[treatment == 0] <- sim_control$event
  }

  sim_treatment <- pwe_sim(
    hazard = hazard_treatment,
    n = sum(treatment == 1),
    maxtime = end_of_study,
    cutpoints = generation_cutpoints
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
