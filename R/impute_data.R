#' @title Impute incomplete data set with event times conditional on observed
#'   data
#'
#' @description Imputes the incomplete time-to-event data with event times from
#'   a piecewise exponential model conditional on what is already observed. This
#'   can be applied to any input data frame, whether for expected success,
#'   futility, or final analysis imputation.
#'
#' @inheritParams survival_adapt
#' @inheritParams sim_comp_data
#' @inheritParams haz_to_prop
#' @param data_in data frame. The time-to-event data that requires imputation,
#'   with columns for treatment assignment (`treatment`, coded `1`
#'   for treatment and `0` for control; single-arm designs use all 1s),
#'   event time (`time`), event indicator (`event`), and indicators
#'   of whether the subject requires imputation for expected success
#'   (`subject_impute_success`) or futility (`subject_impute_futility`).
#' @param hazard array. Hazard parameters for the piecewise exponential
#'   distribution. This should be a single sample from a posterior distribution.
#'   The single slice must have dimensions 1 (rows), \eqn{J} (columns), and 2
#'   (third dimension, in order of treatment and control hazard slices).
#' @param type character. Whether imputation is for `success` or `futility`.
#'
#' @details This function is not intended to be used outside of the main
#'   simulation programs: [survival_adapt()] and [sim_trials()].
#'
#' @return A data frame with the same number of rows (number of subjects), but
#'   with imputed event times.
#'
#' @noRd
impute_data <- function(
  data_in,
  hazard,
  end_of_study,
  cutpoints,
  type,
  single_arm
) {
  stopifnot(dim(hazard)[1] == 1) # otherwise, need to use vectorized ppwe()

  # Start from the original data and update only the rows that need imputation.
  # This keeps the incoming row order and avoids rebuilding the full data frame.
  data_impute <- data_in

  # Pick the imputation flag and simulator for the requested analysis. Success
  # imputations condition on observed follow-up; futility imputations simulate
  # complete outcomes for not-yet-enrolled subjects.
  if (type == "success") {
    subject_requires_imputation <- data_in$subject_impute_success
    impute <- function(idx, hazard_slice) {
      pwe_impute(
        time = data_in$time[idx],
        hazard = hazard[1, , hazard_slice],
        maxtime = end_of_study,
        cutpoints = cutpoints
      )
    }
  } else if (type == "futility") {
    subject_requires_imputation <- data_in$subject_impute_futility
    impute <- function(idx, hazard_slice) {
      pwe_sim(
        n = sum(idx),
        hazard = hazard[1, , hazard_slice],
        maxtime = end_of_study,
        cutpoints = cutpoints
      )
    }
  } else {
    stop("'type' must be either 'success' or 'futility'")
  }

  # Preserve the old RNG order: treatment rows are imputed before control rows.
  # The data column uses treatment = 1 for treatment and treatment = 0 for
  # control; the hazard array uses slice 1 for treatment and slice 2 for control.
  treatment_idx <- data_in$treatment == 1 & subject_requires_imputation
  impute_treatment <- impute(treatment_idx, hazard_slice = 1)
  data_impute$time[treatment_idx] <- impute_treatment$time
  data_impute$event[treatment_idx] <- impute_treatment$event

  # Two-arm studies use the second hazard slice for control-arm imputations.
  if (!single_arm) {
    control_idx <- data_in$treatment == 0 & subject_requires_imputation
    impute_control <- impute(control_idx, hazard_slice = 2)
    data_impute$time[control_idx] <- impute_control$time
    data_impute$event[control_idx] <- impute_control$event
  }

  # Check: imputed data should have same number of subjects as
  #        the interim data
  if (nrow(data_impute) != nrow(data_in)) {
    stop("Number of subjects different after imputation!")
  }

  return(data_impute)
}
