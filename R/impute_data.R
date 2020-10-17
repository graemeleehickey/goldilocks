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
#' @param hazard vector. Hazard parameters for the piecewise exponential
#'   distribution. This will most likely be a single sample from a posterior
#'   distribution.
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
#' @export
impute_data <- function(data_in, hazard, end_of_study, cutpoint, type,
                        single_arm) {

  ## Expected success
  if (type == "success") {

    # Imputing for treatment group
    treatment_impute <- data_in %>%
      filter(treatment == 1 & subject_impute_success)

    impute_treatment <- pwe_impute(time     = treatment_impute$time,
                                   hazard   = hazard,
                                   maxtime  = end_of_study,
                                   cutpoint = cutpoint)

    # Imputing for control group
    if (!single_arm) {
      control_impute <- data_in %>%
        filter(treatment == 0 & subject_impute_success)

      # Impute PWE event times conditional on current observed time
      impute_control <- pwe_impute(time     = control_impute$time,
                                   hazard   = hazard,
                                   maxtime  = end_of_study,
                                   cutpoint = cutpoint)
    }
  ## Futility
  } else if (type == "futility") {

    # Imputing for treatment group
    treatment_impute <- data_in %>%
      filter(treatment == 1 & subject_impute_futility)

    impute_treatment <- pwe_sim(n        = nrow(treatment_impute),
                                hazard   = hazard,
                                maxtime  = end_of_study,
                                cutpoint = cutpoint)

    # Imputing for control group
    if (!single_arm) {
      control_impute <- data_in %>%
        filter(treatment == 0 & subject_impute_futility)

      # Impute PWE event times conditional on current observed time
      impute_control <- pwe_sim(n        = nrow(control_impute),
                                hazard   = hazard,
                                maxtime  = end_of_study,
                                cutpoint = cutpoint)
    }
  }

  data_treatment_impute <- treatment_impute %>%
    bind_cols(time_impute  = impute_treatment$time,
              event_impute = impute_treatment$event)

  if (!single_arm) {
    data_control_impute <- control_impute %>%
      bind_cols(time_impute  = impute_control$time,
                event_impute = impute_control$event)
  } else {
    data_control_impute <- NULL
  }

  # Non-imputed data
  data_noimpute <- data_in %>%
    mutate(time_impute = time,
           event_impute = event)
  if (type == "success") {
    data_noimpute <- data_noimpute %>%
      filter(!subject_impute_success)
  } else if (type == "futility") {
    data_noimpute <- data_noimpute %>%
      filter(!subject_impute_futility)
  }

  # Combine imputed and non-imputed data
  data_impute <- bind_rows(data_control_impute,
                           data_treatment_impute,
                           data_noimpute) %>%
    mutate(time = time_impute,
           event = event_impute) %>%
    select(-c(time_impute, event_impute))

  # Check: imputed data should have same number of subjects as
  #        the interim data
  if (nrow(data_impute) != nrow(data_in)) {
    stop("Number of subjects different after imputation!")
  }

  return(data_impute)

}
