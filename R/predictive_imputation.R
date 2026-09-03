#' Generate a batch of interim predictive imputations
#'
#' @description Generates expected-success and optional futility imputations
#'   for every posterior hazard draw. Only rows that require imputation are
#'   retained; observed outcomes remain unchanged when a completed data set is
#'   constructed for a posterior draw.
#'
#' @param data_in A data frame containing the observed interim data and
#'   indicators identifying outcomes that require predictive imputation.
#' @param hazards A three-dimensional numeric array of posterior hazard draws,
#'   with draws in rows, piecewise intervals in columns, and treatment followed
#'   by control in the third dimension.
#' @inheritParams impute_data
#' @param check_futility A single logical value indicating whether outcomes for
#'   subjects not yet enrolled must be simulated for the futility calculation.
#'
#' @return A list containing subject-by-draw time and event matrices for the
#'   current cohort and, when requested, future subjects.
#'
#' @details Random inputs are generated in posterior-draw order. Within each
#'   draw, current-cohort treatment and control imputations precede future-
#'   subject treatment and control imputations. This matches the scalar
#'   imputation order when completed-data analysis itself consumes no random
#'   numbers. In the batched interim calculation, all imputations precede all
#'   completed-data analyses.
#'
#' @keywords internal
#' @noRd
impute_predictive_draws <- function(
  data_in,
  hazards,
  end_of_study,
  cutpoints,
  single_arm,
  binary_imputation = c("event-time", "bernoulli"),
  check_futility
) {
  binary_imputation <- match.arg(binary_imputation)
  validate_cutpoints(cutpoints)
  validate_endpoint_time(end_of_study, cutpoints, "end_of_study")

  required_columns <- c(
    "time",
    "event",
    "treatment",
    "subject_enrolled",
    "subject_impute_success",
    "subject_impute_futility"
  )
  if (
    !is.data.frame(data_in) ||
      !all(required_columns %in% names(data_in))
  ) {
    stop(
      "Internal imputation invariant failed: incomplete interim data",
      call. = FALSE
    )
  }
  if (
    anyNA(data_in[required_columns]) ||
      any(!data_in$treatment %in% c(0, 1)) ||
      !is.logical(data_in$subject_enrolled) ||
      !is.logical(data_in$subject_impute_success) ||
      !is.logical(data_in$subject_impute_futility)
  ) {
    stop(
      "Internal imputation invariant failed: invalid interim data",
      call. = FALSE
    )
  }
  if (
    length(single_arm) != 1L ||
      !is.logical(single_arm) ||
      is.na(single_arm) ||
      length(check_futility) != 1L ||
      !is.logical(check_futility) ||
      is.na(check_futility)
  ) {
    stop(
      "Internal imputation invariant failed: invalid design flags",
      call. = FALSE
    )
  }
  if (
    any(data_in$subject_impute_success & !data_in$subject_enrolled) ||
      any(data_in$subject_impute_futility & data_in$subject_enrolled) ||
      (single_arm && any(data_in$treatment != 1L))
  ) {
    stop(
      "Internal imputation invariant failed: inconsistent imputation flags",
      call. = FALSE
    )
  }

  hazard_dimensions <- dim(hazards)
  required_arms <- if (single_arm) 1L else 2L
  if (
    !is.array(hazards) ||
      !is.numeric(hazards) ||
      length(hazard_dimensions) != 3L ||
      hazard_dimensions[1L] < 1L ||
      hazard_dimensions[2L] != length(cutpoints) + 1L ||
      hazard_dimensions[3L] < required_arms
  ) {
    stop(
      "Internal imputation invariant failed: invalid posterior hazard array",
      call. = FALSE
    )
  }
  relevant_hazards <- hazards[,,
    seq_len(required_arms),
    drop = FALSE
  ]
  if (
    anyNA(relevant_hazards) ||
      any(!is.finite(relevant_hazards)) ||
      any(relevant_hazards < 0)
  ) {
    stop(
      "Internal imputation invariant failed: invalid posterior hazards",
      call. = FALSE
    )
  }

  n_draws <- hazard_dimensions[1L]
  current_flag <- data_in$subject_impute_success
  future_flag <- data_in$subject_impute_futility & check_futility
  current <- empty_predictive_imputations(which(current_flag), n_draws)
  future <- if (check_futility) {
    empty_predictive_imputations(which(future_flag), n_draws)
  } else {
    NULL
  }

  current_treatment <- which(current_flag & data_in$treatment == 1L)
  current_control <- if (single_arm) {
    integer()
  } else {
    which(current_flag & data_in$treatment == 0L)
  }
  future_treatment <- which(future_flag & data_in$treatment == 1L)
  future_control <- if (single_arm) {
    integer()
  } else {
    which(future_flag & data_in$treatment == 0L)
  }

  random_inputs <- list(
    current_treatment = matrix(
      numeric(length(current_treatment) * n_draws),
      nrow = length(current_treatment),
      ncol = n_draws
    ),
    current_control = matrix(
      numeric(length(current_control) * n_draws),
      nrow = length(current_control),
      ncol = n_draws
    ),
    future_treatment = matrix(
      numeric(length(future_treatment) * n_draws),
      nrow = length(future_treatment),
      ncol = n_draws
    ),
    future_control = matrix(
      numeric(length(future_control) * n_draws),
      nrow = length(future_control),
      ncol = n_draws
    )
  )

  for (draw in seq_len(n_draws)) {
    random_inputs$current_treatment[, draw] <-
      predictive_random_values(
        length(current_treatment),
        type = "success",
        binary_imputation = binary_imputation,
        cutpoints = cutpoints
      )
    random_inputs$current_control[, draw] <-
      predictive_random_values(
        length(current_control),
        type = "success",
        binary_imputation = binary_imputation,
        cutpoints = cutpoints
      )
    random_inputs$future_treatment[, draw] <-
      predictive_random_values(
        length(future_treatment),
        type = "futility",
        binary_imputation = binary_imputation,
        cutpoints = cutpoints
      )
    random_inputs$future_control[, draw] <-
      predictive_random_values(
        length(future_control),
        type = "futility",
        binary_imputation = binary_imputation,
        cutpoints = cutpoints
      )
  }

  treatment_hazards <- matrix(
    hazards[,, 1L],
    nrow = n_draws,
    ncol = hazard_dimensions[2L]
  )
  current <- fill_predictive_imputations(
    imputations = current,
    rows = current_treatment,
    values = impute_predictive_arm(
      time = data_in$time[current_treatment],
      hazards = treatment_hazards,
      random_input = random_inputs$current_treatment,
      end_of_study = end_of_study,
      cutpoints = cutpoints,
      binary_imputation = binary_imputation
    )
  )
  if (check_futility) {
    future <- fill_predictive_imputations(
      imputations = future,
      rows = future_treatment,
      values = impute_predictive_arm(
        time = rep.int(0, length(future_treatment)),
        hazards = treatment_hazards,
        random_input = random_inputs$future_treatment,
        end_of_study = end_of_study,
        cutpoints = cutpoints,
        binary_imputation = binary_imputation
      )
    )
  }

  if (!single_arm) {
    control_hazards <- matrix(
      hazards[,, 2L],
      nrow = n_draws,
      ncol = hazard_dimensions[2L]
    )
    current <- fill_predictive_imputations(
      imputations = current,
      rows = current_control,
      values = impute_predictive_arm(
        time = data_in$time[current_control],
        hazards = control_hazards,
        random_input = random_inputs$current_control,
        end_of_study = end_of_study,
        cutpoints = cutpoints,
        binary_imputation = binary_imputation
      )
    )
    if (check_futility) {
      future <- fill_predictive_imputations(
        imputations = future,
        rows = future_control,
        values = impute_predictive_arm(
          time = rep.int(0, length(future_control)),
          hazards = control_hazards,
          random_input = random_inputs$future_control,
          end_of_study = end_of_study,
          cutpoints = cutpoints,
          binary_imputation = binary_imputation
        )
      )
    }
  }

  list(
    n_draws = n_draws,
    current = current,
    future = future,
    rng_order = paste(
      "posterior draw; current treatment, current control,",
      "future treatment, future control"
    )
  )
}

#' Create empty predictive imputations
#'
#' @keywords internal
#' @noRd
empty_predictive_imputations <- function(rows, n_draws) {
  list(
    rows = as.integer(rows),
    time = matrix(
      numeric(length(rows) * n_draws),
      nrow = length(rows),
      ncol = n_draws
    ),
    event = matrix(
      integer(length(rows) * n_draws),
      nrow = length(rows),
      ncol = n_draws
    )
  )
}

#' Draw the random values used for predictive imputation
#'
#' @keywords internal
#' @noRd
predictive_random_values <- function(
  n,
  type,
  binary_imputation,
  cutpoints
) {
  if (n == 0L) {
    return(numeric())
  }
  if (binary_imputation == "bernoulli") {
    return(runif(n))
  }
  if (type == "success") {
    return(-log1p(-runif(n)))
  }
  if (length(cutpoints) == 0L) {
    return(rexp(n, rate = 1))
  }
  -log(runif(n))
}

#' Impute one arm across posterior hazard draws
#'
#' @keywords internal
#' @noRd
impute_predictive_arm <- function(
  time,
  hazards,
  random_input,
  end_of_study,
  cutpoints,
  binary_imputation
) {
  n_subjects <- length(time)
  n_draws <- nrow(hazards)
  if (n_subjects == 0L) {
    return(list(
      time = matrix(numeric(), nrow = 0L, ncol = n_draws),
      event = matrix(integer(), nrow = 0L, ncol = n_draws)
    ))
  }

  interval_lower <- c(0, cutpoints)
  interval_upper <- c(cutpoints, Inf)
  if (binary_imputation == "bernoulli") {
    remaining_exposure <- outer(
      time,
      seq_len(ncol(hazards)),
      function(observed_time, interval) {
        pmax(
          0,
          pmin(end_of_study, interval_upper[interval]) -
            pmax(observed_time, interval_lower[interval])
        )
      }
    )
    probability <- cumulative_hazard_to_probability(
      remaining_exposure %*% t(hazards)
    )
    return(list(
      time = matrix(
        rep.int(end_of_study, n_subjects * n_draws),
        nrow = n_subjects,
        ncol = n_draws
      ),
      event = matrix(
        as.integer(random_input <= probability),
        nrow = n_subjects,
        ncol = n_draws
      )
    ))
  }

  exposure_to_time <- outer(
    time,
    seq_len(ncol(hazards)),
    function(observed_time, interval) {
      pmax(
        0,
        pmin(observed_time, interval_upper[interval]) -
          interval_lower[interval]
      )
    }
  )
  target_hazard <- exposure_to_time %*% t(hazards) + random_input
  interval_widths <- endpoint_interval_widths(cutpoints, end_of_study)
  cumulative_hazard <- sweep(hazards, 2L, interval_widths, `*`)
  if (ncol(cumulative_hazard) > 1L) {
    for (interval in 2:ncol(cumulative_hazard)) {
      cumulative_hazard[, interval] <-
        cumulative_hazard[, interval - 1L] +
        cumulative_hazard[, interval]
    }
  }

  imputed_time <- matrix(
    rep.int(end_of_study, n_subjects * n_draws),
    nrow = n_subjects,
    ncol = n_draws
  )
  imputed_event <- matrix(
    0L,
    nrow = n_subjects,
    ncol = n_draws
  )
  unresolved <- matrix(
    TRUE,
    nrow = n_subjects,
    ncol = n_draws
  )
  observed_time <- matrix(
    time,
    nrow = n_subjects,
    ncol = n_draws
  )

  for (interval in seq_len(ncol(hazards))) {
    cumulative_before <- if (interval == 1L) {
      numeric(n_draws)
    } else {
      cumulative_hazard[, interval - 1L]
    }
    cumulative_end <- matrix(
      cumulative_hazard[, interval],
      nrow = n_subjects,
      ncol = n_draws,
      byrow = TRUE
    )
    cumulative_before <- matrix(
      cumulative_before,
      nrow = n_subjects,
      ncol = n_draws,
      byrow = TRUE
    )
    interval_hazard <- matrix(
      hazards[, interval],
      nrow = n_subjects,
      ncol = n_draws,
      byrow = TRUE
    )
    selected <- unresolved &
      interval_hazard > 0 &
      target_hazard <= cumulative_end
    if (any(selected)) {
      candidate_time <- interval_lower[interval] +
        (target_hazard - cumulative_before) / interval_hazard
      candidate_time <- pmin(
        end_of_study,
        pmax(observed_time, candidate_time)
      )
      imputed_time[selected] <- candidate_time[selected]
      imputed_event[selected] <- 1L
      unresolved[selected] <- FALSE
    }
  }

  list(time = imputed_time, event = imputed_event)
}

#' Add one arm to the predictive imputations
#'
#' @keywords internal
#' @noRd
fill_predictive_imputations <- function(imputations, rows, values) {
  if (length(rows) == 0L) {
    return(imputations)
  }
  positions <- match(rows, imputations$rows)
  imputations$time[positions, ] <- values$time
  imputations$event[positions, ] <- values$event
  imputations
}
