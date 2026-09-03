#' Complete one predictive dataset
#'
#' @keywords internal
#' @noRd
complete_predictive_data <- function(
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
prepare_predictive_outcomes <- function(
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

#' Complete one predictive survival outcome
#'
#' @description Inserts one predictive replicate into the outcome vectors
#'   prepared by `prepare_predictive_outcomes()`. Treatment assignments
#'   are reused because they do not vary across predictive replicates.
#'
#' @param prepared Prepared predictive outcomes returned by
#'   `prepare_predictive_outcomes()`.
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
complete_predictive_outcomes <- function(
  prepared,
  imputations,
  draw,
  maximum = FALSE
) {
  base <- if (maximum) prepared$maximum else prepared$current
  if (is.null(base)) {
    stop(
      paste0(
        "Internal predictive-analysis invariant failed: maximum-sample ",
        "outcomes unavailable"
      ),
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
predictive_binary_counts <- function(
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
