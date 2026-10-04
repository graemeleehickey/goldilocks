# Calendar-masked data for an enrollment-triggered interim analysis.
prepare_simulated_interim <- function(data_total, planned_N, end_of_study) {
  data_interim <- within(data_total, {
    subject_enrolled <- (id <= planned_N)
    subject_impute_futility <- !subject_enrolled
    time_from_rand_at_look <- enrollment[planned_N] -
      enrollment
    subject_impute_success <-
      # Had event, but has not occurred yet (based on interim look)
      ((event == 1) * (time_from_rand_at_look < time) & subject_enrolled) |
      # Event-free and not had opportunity to complete full follow
      ((event == 0) *
        (time_from_rand_at_look < end_of_study) &
        subject_enrolled) |
      (loss_to_fu & subject_enrolled)
  })

  # Mask the data at time of look
  # Note: subjects at the exact interim boundary have
  # time_from_rand_at_look = 0, yielding time = 0 after masking.
  # Clamp to .Machine$double.eps so the boundary subject contributes
  # negligible but non-zero exposure to the interim posterior.
  data_interim <- within(data_interim, {
    time <- pmax(pmin(time, time_from_rand_at_look), .Machine$double.eps)
    event <- ifelse(subject_impute_success, 0, event)
  })

  data_interim
}
