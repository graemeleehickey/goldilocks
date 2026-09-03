#' @title Complete one data set from a posterior hazard draw
#'
#' @description Imputes incomplete time-to-event outcomes conditional on the
#'   observed follow-up and one draw from the posterior hazard distribution.
#'   This function remains the active calculation for final-analysis imputation
#'   and provides the single-imputation reference calculation used in tests.
#'   Interim analyses generate all predictive imputations together using
#'   `impute_predictive_draws()`.
#'
#' @inheritParams survival_adapt
#' @inheritParams sim_comp_data
#' @inheritParams haz_to_prop
#' @param data_in A data frame with one row per subject and columns for
#'   treatment assignment (`treatment`, coded `1`
#'   for treatment and `0` for control; single-arm designs use all 1s),
#'   event time (`time`), event indicator (`event`), and indicators
#'   of whether the subject requires imputation for expected success
#'   (`subject_impute_success`) or futility (`subject_impute_futility`).
#' @param hazard A three-dimensional numeric array containing one posterior draw
#'   of the piecewise-exponential hazard rates. Its dimensions must be `1` by
#'   \eqn{J} by `2`, where \eqn{J} is the number of intervals and the third
#'   dimension contains treatment followed by control. The control values are
#'   ignored for a single-arm design.
#' @param type A required character string: `"success"` completes outcomes for
#'   currently enrolled subjects, whereas `"futility"` simulates outcomes for
#'   subjects who could be enrolled up to the maximum sample size.
#' @param binary_imputation A character string specifying whether to impute a
#'   conditional event time (`"event-time"`, the default) or draw the endpoint
#'   event indicator directly (`"bernoulli"`).
#'
#' @details This is an active internal function, not an archived function. It is
#'   used by the final-analysis calculations in [survival_adapt()] and
#'   [sim_trials()].
#'
#' @return A data frame with the same rows and columns as `data_in`, with the
#'   required `time` and `event` values replaced by imputed outcomes.
#'
#' @noRd
impute_data <- function(
  data_in,
  hazard,
  end_of_study,
  cutpoints,
  type,
  single_arm,
  binary_imputation = c("event-time", "bernoulli")
) {
  binary_imputation <- match.arg(binary_imputation)

  if (
    !is.array(hazard) ||
      length(dim(hazard)) != 3 ||
      dim(hazard)[1] != 1
  ) {
    stop(
      "'hazard' must be a three-dimensional array with exactly one posterior draw"
    )
  }

  # Start from the original data and update only the rows that need imputation.
  # This keeps the incoming row order and avoids rebuilding the full data frame.
  data_impute <- data_in

  # Pick the imputation flag. Success imputations condition on observed
  # follow-up; futility imputations simulate complete outcomes for
  # not-yet-enrolled subjects.
  if (type == "success") {
    subject_requires_imputation <- data_in$subject_impute_success
  } else if (type == "futility") {
    subject_requires_imputation <- data_in$subject_impute_futility
  } else {
    stop("'type' must be either 'success' or 'futility'")
  }

  if (binary_imputation == "bernoulli") {
    impute <- function(idx, hazard_slice) {
      conditioning_time <- if (type == "success") {
        data_in$time[idx]
      } else {
        rep.int(0, sum(idx))
      }
      event_probability <- pwe_conditional_event_probability(
        time = conditioning_time,
        hazard = hazard[1, , hazard_slice],
        end_of_study = end_of_study,
        cutpoints = cutpoints
      )
      data.frame(
        # The precise event time is not drawn in this mode. Binary analyses use
        # only endpoint status, so end_of_study represents completed follow-up.
        time = rep.int(end_of_study, length(event_probability)),
        event = as.numeric(
          runif(length(event_probability)) <= event_probability
        )
      )
    }
  } else if (type == "success") {
    impute <- function(idx, hazard_slice) {
      pwe_impute(
        time = data_in$time[idx],
        hazard = hazard[1, , hazard_slice],
        maxtime = end_of_study,
        cutpoints = cutpoints
      )
    }
  } else {
    impute <- function(idx, hazard_slice) {
      pwe_sim(
        n = sum(idx),
        hazard = hazard[1, , hazard_slice],
        maxtime = end_of_study,
        cutpoints = cutpoints
      )
    }
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
