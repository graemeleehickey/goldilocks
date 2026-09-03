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
  current <- new_predictive_imputation_block(which(current_flag), n_draws)
  future <- if (check_futility) {
    new_predictive_imputation_block(which(future_flag), n_draws)
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
      predictive_random_input(
        length(current_treatment),
        type = "success",
        binary_imputation = binary_imputation,
        cutpoints = cutpoints
      )
    random_inputs$current_control[, draw] <-
      predictive_random_input(
        length(current_control),
        type = "success",
        binary_imputation = binary_imputation,
        cutpoints = cutpoints
      )
    random_inputs$future_treatment[, draw] <-
      predictive_random_input(
        length(future_treatment),
        type = "futility",
        binary_imputation = binary_imputation,
        cutpoints = cutpoints
      )
    random_inputs$future_control[, draw] <-
      predictive_random_input(
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
  current <- fill_predictive_imputation_block(
    block = current,
    rows = current_treatment,
    values = impute_predictive_group(
      time = data_in$time[current_treatment],
      hazards = treatment_hazards,
      random_input = random_inputs$current_treatment,
      end_of_study = end_of_study,
      cutpoints = cutpoints,
      binary_imputation = binary_imputation
    )
  )
  if (check_futility) {
    future <- fill_predictive_imputation_block(
      block = future,
      rows = future_treatment,
      values = impute_predictive_group(
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
    current <- fill_predictive_imputation_block(
      block = current,
      rows = current_control,
      values = impute_predictive_group(
        time = data_in$time[current_control],
        hazards = control_hazards,
        random_input = random_inputs$current_control,
        end_of_study = end_of_study,
        cutpoints = cutpoints,
        binary_imputation = binary_imputation
      )
    )
    if (check_futility) {
      future <- fill_predictive_imputation_block(
        block = future,
        rows = future_control,
        values = impute_predictive_group(
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

#' Allocate storage for a predictive-imputation block
#'
#' @keywords internal
#' @noRd
new_predictive_imputation_block <- function(rows, n_draws) {
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

#' Draw the scalar random input used by one imputation group
#'
#' @keywords internal
#' @noRd
predictive_random_input <- function(
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
impute_predictive_group <- function(
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

#' Fill one arm in a predictive-imputation block
#'
#' @keywords internal
#' @noRd
fill_predictive_imputation_block <- function(block, rows, values) {
  if (length(rows) == 0L) {
    return(block)
  }
  positions <- match(rows, block$rows)
  block$time[positions, ] <- values$time
  block$event[positions, ] <- values$event
  block
}

#' Materialize one completed predictive dataset
#'
#' @keywords internal
#' @noRd
materialize_predictive_draw <- function(
  data_in,
  imputations,
  draw,
  include_future = FALSE
) {
  if (
    length(draw) != 1L ||
      !is.numeric(draw) ||
      is.na(draw) ||
      draw != as.integer(draw) ||
      draw < 1L ||
      draw > imputations$n_draws
  ) {
    stop("Internal imputation invariant failed: invalid draw index")
  }
  if (include_future && is.null(imputations$future)) {
    stop("Internal imputation invariant failed: future draws unavailable")
  }

  data_out <- data_in
  current_rows <- imputations$current$rows
  data_out$time[current_rows] <- imputations$current$time[, draw]
  data_out$event[current_rows] <- imputations$current$event[, draw]
  if (include_future) {
    future_rows <- imputations$future$rows
    data_out$time[future_rows] <- imputations$future$time[, draw]
    data_out$event[future_rows] <- imputations$future$event[, draw]
  }
  data_out
}

#' Prepare outcome vectors for repeated survival analyses
#'
#' @description Retains the observed `time`, `event`, and `treatment` vectors
#'   needed for completed-data survival analyses and calculates the positions
#'   at which current and future imputations will be inserted. This preparation
#'   is performed once for an interim look.
#'
#' @param data_in A data frame containing the prepared interim data.
#' @param imputations Predictive imputations returned by
#'   `impute_predictive_draws()`.
#' @param check_futility A single logical value indicating whether
#'   maximum-sample-size outcomes are required.
#'
#' @return A list containing fixed current- and maximum-sample outcome vectors
#'   and the corresponding imputation positions.
#'
#' @keywords internal
#' @noRd
prepare_predictive_survival_state <- function(
  data_in,
  imputations,
  check_futility
) {
  required_columns <- c(
    "time",
    "event",
    "treatment",
    "subject_enrolled"
  )
  missing_columns <- setdiff(required_columns, names(data_in))
  if (!is.data.frame(data_in) || length(missing_columns) > 0L) {
    stop(
      "Internal predictive-analysis invariant failed: incomplete interim data",
      call. = FALSE
    )
  }
  if (
    !is.list(imputations) ||
      !is.numeric(imputations$n_draws) ||
      length(imputations$n_draws) != 1L ||
      imputations$n_draws < 1L ||
      is.null(imputations$current)
  ) {
    stop(
      "Internal predictive-analysis invariant failed: invalid imputations",
      call. = FALSE
    )
  }
  if (
    length(check_futility) != 1L ||
      !is.logical(check_futility) ||
      is.na(check_futility) ||
      (check_futility && is.null(imputations$future))
  ) {
    stop(
      "Internal predictive-analysis invariant failed: invalid futility state",
      call. = FALSE
    )
  }

  enrolled_rows <- which(data_in$subject_enrolled)
  current_positions <- match(imputations$current$rows, enrolled_rows)
  if (anyNA(current_positions)) {
    stop(
      paste0(
        "Internal predictive-analysis invariant failed: current imputation ",
        "outside the enrolled cohort"
      ),
      call. = FALSE
    )
  }

  current <- list(
    time = data_in$time[enrolled_rows],
    event = data_in$event[enrolled_rows],
    treatment = data_in$treatment[enrolled_rows],
    current_positions = current_positions,
    future_positions = integer()
  )
  maximum <- if (check_futility) {
    list(
      time = data_in$time,
      event = data_in$event,
      treatment = data_in$treatment,
      current_positions = imputations$current$rows,
      future_positions = imputations$future$rows
    )
  } else {
    NULL
  }

  list(current = current, maximum = maximum)
}

#' Complete one predictive survival outcome from stored vectors
#'
#' @description Inserts one predictive replicate into the outcome vectors
#'   prepared by `prepare_predictive_survival_state()`. Treatment assignments
#'   are reused because they do not vary across predictive replicates.
#'
#' @param state A prepared predictive survival state.
#' @param imputations Predictive imputations returned by
#'   `impute_predictive_draws()`.
#' @param draw A positive integer identifying the predictive replicate.
#' @param maximum A single logical value indicating whether to include future
#'   participants and analyze the maximum sample size.
#'
#' @return A list containing completed `time`, `event`, and `treatment` vectors.
#'
#' @keywords internal
#' @noRd
complete_predictive_survival_draw <- function(
  state,
  imputations,
  draw,
  maximum = FALSE
) {
  base <- if (maximum) state$maximum else state$current
  if (is.null(base)) {
    stop(
      "Internal predictive-analysis invariant failed: maximum state unavailable",
      call. = FALSE
    )
  }

  time <- base$time
  event <- base$event
  time[base$current_positions] <- imputations$current$time[, draw]
  event[base$current_positions] <- imputations$current$event[, draw]
  if (maximum) {
    time[base$future_positions] <- imputations$future$time[, draw]
    event[base$future_positions] <- imputations$future$event[, draw]
  }

  list(
    time = time,
    event = event,
    treatment = base$treatment
  )
}

#' Summarize completed binary outcomes for predictive replicates
#'
#' @description Combines observed events with the imputed endpoint statuses of
#'   pending and future subjects. For each predictive replicate, it records the
#'   number of events and subjects in each arm. These are the sufficient
#'   statistics for the risk-difference and beta-binomial analyses.
#'
#' @param data_in A data frame containing the observed event indicators,
#'   treatment assignments, enrollment status, and indicators of which current
#'   or future outcomes require imputation.
#' @param imputations A list of current-cohort and future-subject endpoint
#'   imputations returned by `impute_predictive_draws()`.
#' @param single_arm A single logical value indicating a one-arm design.
#' @param check_futility A single logical value indicating whether maximum-
#'   sample-size completions are required for the futility calculation.
#'
#' @return A list containing data frames for the currently enrolled cohort and,
#'   when requested, the maximum sample size. Each row represents one
#'   predictive replicate and contains control and treatment event counts and
#'   sample sizes.
#'
#' @keywords internal
#' @noRd
predictive_binary_count_states <- function(
  data_in,
  imputations,
  single_arm,
  check_futility
) {
  n_draws <- imputations$n_draws
  arm_values <- if (single_arm) 1L else c(0L, 1L)

  block_event_counts <- function(block) {
    counts <- matrix(
      0L,
      nrow = length(arm_values),
      ncol = n_draws,
      dimnames = list(as.character(arm_values), NULL)
    )
    if (length(block$rows) == 0L) {
      return(counts)
    }
    block_treatment <- data_in$treatment[block$rows]
    for (arm in arm_values) {
      rows <- block_treatment == arm
      if (any(rows)) {
        counts[as.character(arm), ] <- colSums(
          block$event[rows, , drop = FALSE]
        )
      }
    }
    counts
  }

  fixed_rows <- data_in$subject_enrolled &
    !data_in$subject_impute_success
  fixed_events <- vapply(
    arm_values,
    function(arm) {
      sum(data_in$event[fixed_rows & data_in$treatment == arm])
    },
    numeric(1L)
  )
  names(fixed_events) <- as.character(arm_values)

  current_events <- sweep(
    block_event_counts(imputations$current),
    MARGIN = 1L,
    STATS = fixed_events,
    FUN = `+`
  )
  count_frame <- function(events, maximum = FALSE) {
    n_for_arm <- function(arm) {
      rows <- data_in$treatment == arm
      if (!maximum) {
        rows <- rows & data_in$subject_enrolled
      }
      sum(rows)
    }
    data.frame(
      events_control = if (single_arm) {
        rep.int(0L, n_draws)
      } else {
        as.integer(events["0", ])
      },
      n_control = if (single_arm) {
        rep.int(0L, n_draws)
      } else {
        rep.int(as.integer(n_for_arm(0L)), n_draws)
      },
      events_treatment = as.integer(events["1", ]),
      n_treatment = rep.int(as.integer(n_for_arm(1L)), n_draws),
      stringsAsFactors = FALSE
    )
  }

  current <- count_frame(current_events)
  maximum <- if (check_futility) {
    if (is.null(imputations$future)) {
      stop("Internal imputation invariant failed: future draws unavailable")
    }
    count_frame(
      current_events + block_event_counts(imputations$future),
      maximum = TRUE
    )
  } else {
    NULL
  }

  list(current = current, maximum = maximum)
}
