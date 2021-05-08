#' @title Impute incomplete data set with event times conditional on observed
#'   data
#'
#' @description Imputes the incomplete time-to-event data with event times from
#'   a piecewise exponential model conditional on what is already observed. This
#'   can be applied to any input data frame, whether for expected success,
#'   futility, or final analysis imputation.
#'
#' @inheritParams survival_adapt
#' @inheritParams haz_to_prop
#' @param data_in data frame. The time-to-event data that requires imputation,
#'   with columns for treatment arm (\code{treatment}, which can be all 1s if
#'   single-arm), event time (\code{time}), event indicator (\code{event}),
#'   indicator of whether the subject requires imputation for expected success
#'   (\code{subject_impute_success}) or futility
#'   (\code{subject_impute_futility}).
#' @param hazard array. Hazard parameters for the piecewise exponential
#'   distribution. This should be a single sample from a posterior distribution.
#'   The single slice must have dimensions 1 (rows), \eqn{J} (columns), and 2
#'   (third dimension, in order of treatment and control).
#' @param type character. Whether imputation is for \code{success} or
#'   \code{futility}.
#'
#' @details This function is not intended to be used outside of the main
#'   simulation programs: \code{\link{survival_adapt}} and
#'   \code{\link{sim_trials}}.
#'
#' @return A data frame with the same number of rows (number of subjects), but
#'   with imputed event times.
#'
#' @noRd
impute_data <- function(data_in, hazard, end_of_study, cutpoints, type,
                        single_arm) {

  stopifnot(dim(hazard)[1] == 1) # otherwise, need to use vectorized ppwe()

  ## Expected success
  if (type == "success") {

    # Imputing for treatment group
    treatment_impute <- subset(data_in,
                               treatment == 1 & subject_impute_success)

    impute_treatment <- pwe_impute(time      = treatment_impute$time,
                                   hazard    = hazard[1, , 1],
                                   maxtime   = end_of_study,
                                   cutpoints = cutpoints)

    # Imputing for control group
    if (!single_arm) {
      control_impute <- subset(data_in,
                               treatment == 0 & subject_impute_success)

      # Impute PWE event times conditional on current observed time
      impute_control <- pwe_impute(time      = control_impute$time,
                                   hazard    = hazard[1, , 2],
                                   maxtime   = end_of_study,
                                   cutpoints = cutpoints)
    }
    ## Futility
  } else if (type == "futility") {

    # Imputing for treatment group
    treatment_impute <- subset(data_in,
                               treatment == 1 & subject_impute_futility)

    impute_treatment <- pwe_sim(n         = nrow(treatment_impute),
                                hazard    = hazard[1, , 1],
                                maxtime   = end_of_study,
                                cutpoints = cutpoints)

    # Imputing for control group
    if (!single_arm) {
      control_impute <- subset(data_in,
                               treatment == 0 & subject_impute_futility)

      # Impute PWE event times conditional on current observed time
      impute_control <- pwe_sim(n         = nrow(control_impute),
                                hazard    = hazard[1, , 2],
                                maxtime   = end_of_study,
                                cutpoints = cutpoints)
    }
  }

  data_treatment_impute <- cbind(treatment_impute,
                                 "time_impute"  = impute_treatment$time,
                                 "event_impute" = impute_treatment$event)

  if (!single_arm) {
    data_control_impute <- cbind(control_impute,
                                 "time_impute"  = impute_control$time,
                                 "event_impute" = impute_control$event)
  } else {
    data_control_impute <- NULL
  }

  # Non-imputed data
  data_noimpute              <- data_in
  data_noimpute$time_impute  <- data_noimpute$time
  data_noimpute$event_impute <- data_noimpute$event
  if (type == "success") {
    data_noimpute <- subset(data_noimpute, !subject_impute_success)
  } else if (type == "futility") {
    data_noimpute <- subset(data_noimpute, !subject_impute_futility)
  }

  # Combine imputed and non-imputed data
  data_impute <- rbind(data_control_impute,
                       data_treatment_impute,
                       data_noimpute)
  data_impute$time  <- data_impute$time_impute
  data_impute$event <- data_impute$event_impute
  data_impute <- data_impute[, 1:10]

  # Check: imputed data should have same number of subjects as
  #        the interim data
  if (nrow(data_impute) != nrow(data_in)) {
    stop("Number of subjects different after imputation!")
  }

  return(data_impute)

}
