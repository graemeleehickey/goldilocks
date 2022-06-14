#' @title Simulate a complete clinical trial with event data drawn from a
#'   piecewise exponential distribution
#'
#' @param hazard_treatment vector. Constant hazard rates under the treatment
#'   arm.
#' @param hazard_control vector. Constant hazard rates under the control arm.
#' @param cutpoints  vector. Times at which the baseline hazard changes. Default
#'   is \code{cutpoints = 0}, which corresponds to a simple (non-piecewise)
#'   exponential model.
#' @param N_total integer. Maximum sample size allowable
#' @param lambda vector. Enrollment rates across simulated enrollment times. See
#'   \code{\link{enrollment}} for more details.
#' @param lambda_time vector. Enrollment time(s) at which the enrollment rates
#'   change. Must be same length as lambda. See \code{\link{enrollment}} for
#'   more details.
#' @param end_of_study scalar. Length of the study; i.e. time at which endpoint
#'   will be evaluated.
#' @param block scalar. Block size for generating the randomization schedule.
#' @param rand_ratio vector. Randomization allocation for the ratio of control
#'   to treatment. Integer values mapping the size of the block. See
#'   \code{\link{randomization}} for more details.
#' @param prop_loss scalar. Overall proportion of subjects lost to follow-up.
#'   Defaults to zero.
#'
#' @return A data frame with 1 row per subject and columns:
#' \describe{
#'     \item{\code{time:}}{
#'       numeric. Time of event or censoring time.
#'     }
#'     \item{\code{treatment:}}{
#'       integer. Treatment arm with values \code{1L} for experimental arm, and
#'       \code{0L} for control arm (only if \code{hazard_control} is given).
#'     }
#'     \item{\code{event:}}{
#'       integer. Indicator of whether event occurred (\code{=1L} if occurred
#'       and \code{=0L} if right-censored).
#'     }
#'     \item{\code{enrollment:}}{
#'       numeric. Time of patient enrollment relative to time trial enrolled
#'       first patient.
#'     }
#'     \item{\code{id:}}{
#'       integer. Identification number for each patient.
#'     }
#'     \item{\code{loss_to_fu:}}{
#'       logical. Indicator of whether the patient was lost to follow-up during
#'       the course of observation.
#'     }
#' }
#'
#' @importFrom stats runif sd
#' @export
sim_comp_data <- function(
    hazard_treatment,
    hazard_control    = NULL,
    cutpoints         = 0,
    N_total,
    lambda            = 0.3,
    lambda_time       = 0,
    end_of_study,
    block             = 2,
    rand_ratio        = c(1, 1),
    prop_loss         = 0
    ) {

  ##############################################################################
  ### Run checks on arguments
  ##############################################################################

  # Check: none of the 'cutpoints' exceed 'end_of_study'
  if (!is.null(cutpoints)) {
    stopifnot(any(cutpoints < end_of_study))
  }

  # Assign: indicator of whether single-arm study
  single_arm <- is.null(hazard_control)

  ##############################################################################
  ### Simulate enrollment + assignment
  ##############################################################################

  # Simulate enrollment times
  enrollment <- enrollment(lambda      = lambda,
                           N_total     = N_total,
                           lambda_time = lambda_time)
  enrollment <- enrollment + runif(length(enrollment))
  enrollment <- sort(enrollment)

  # Simulate treatment arm assignment
  if (!single_arm) {
    group <- randomization(N_total    = N_total,
                           block      = block,
                           allocation = rand_ratio)
  } else {
    group <- rep(1, N_total)
  }

  ##############################################################################
  ### Simulate event times
  ##############################################################################

  time  <- rep(NA, length = N_total)
  event <- rep(NA, length = N_total)

  # Simulate TTE outcome
  # - Note: time = time *from* enrollment
  if (!single_arm) {
    sim_control <- pwe_sim(hazard     = hazard_control,
                           n          = sum(!group),
                           maxtime    = end_of_study,
                           cutpoints  = cutpoints)
    time[group == 0]  <- sim_control$time
    event[group == 0] <- sim_control$event
  }

  sim_treatment <- pwe_sim(hazard    = hazard_treatment,
                           n         = sum(group),
                           maxtime   = end_of_study,
                           cutpoints = cutpoints)
  time[group == 1]  <- sim_treatment$time
  event[group == 1] <- sim_treatment$event

  # Simulate loss to follow-up
  loss_to_fu <- rep(FALSE, N_total)
  if (prop_loss > 0) {
    n_loss_to_fu <- ceiling(prop_loss * N_total)
    loss_to_fu[sample(1:N_total, n_loss_to_fu)] <- TRUE
  }

  # Creating a new data.frame for all the variables
  data_total <- data.frame(
    time       = time,
    treatment  = group,
    event      = event,
    enrollment = enrollment,
    id         = 1:N_total,
    loss_to_fu = loss_to_fu)

  # Subjects lost are uniformly distributed
  if (prop_loss > 0) {
    data_total$time[data_total$loss_to_fu]  <- runif(
      n_loss_to_fu, 0, data_total$time[data_total$loss_to_fu])
    data_total$event[data_total$loss_to_fu] <- rep(0, n_loss_to_fu)
  }

  return(data_total)

}
