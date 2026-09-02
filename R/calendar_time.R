#' Summarize operating characteristics on the calendar-time scale
#'
#' @description
#' Summarizes trial duration, follow-up burden, and (when retained) interim-look
#' timing without adding any new simulation arguments. Calendar time is measured
#' from the first patient's enrollment at time zero. `analysis_ready_time` is
#' the time at which the last enrolled subject's observed event or censoring
#' becomes available; it does not include external data-cleaning or database-lock
#' delays.
#'
#' Pass a complete result from [sim_trials()] to retain requested, analyzed, and
#' failed simulation counts. Interim timing requires simulations run with
#' `return_trace = TRUE`. A simulation data frame can also be supplied, but only
#' the trial-duration table can then be calculated.
#'
#' @param data A required complete result returned by [sim_trials()], a
#'   simulation
#'   `data.frame`, or a list of either form. Named list elements identify
#'   scenarios. Existing `scenario` columns and grouping variables are
#'   preserved.
#'
#' @return An object of class `goldilocks_calendar_summary`, containing two wide
#'   data frames:
#'
#'   - `trial_duration`: one row per scenario and stopping reason, plus an
#'     overall row. It reports simulation denominators, stopping counts, sample
#'     size, accrual-stop time, analysis-ready time, planned completion time,
#'     total person-time under follow-up, and peak concurrent follow-up.
#'   - `interim_timing`: one row per scenario and interim look. It reports how
#'     often the look was reached, its calendar time, and the number of subjects
#'     actively under follow-up. This table has zero rows when traces are not
#'     available.
#'
#'   Continuous quantities are reported as means, Monte Carlo standard errors,
#'   and the 10th, 50th, and 90th percentiles. Percentages use the requested
#'   number of simulations as their denominator, so excluded failures remain
#'   visible.
#'
#' @export
summarise_calendar_time <- function(data) {
  normalized <- normalize_summarise_sims_input(data)
  sims <- normalized$data
  required <- c(
    "N_enrolled",
    "accrual_stop_time",
    "analysis_ready_time",
    "planned_completion_time",
    "followup_person_time",
    "peak_active_followup"
  )
  missing <- setdiff(required, names(sims))
  if (length(missing) > 0L) {
    stop(
      "'data' is missing required calendar-time column(s): ",
      paste(missing, collapse = ", "),
      ". Rerun the simulations with this version of goldilocks.",
      call. = FALSE
    )
  }

  if (!"stopping_reason" %in% names(sims)) {
    stopping_columns <- c("stop_futility", "stop_expected_success")
    missing_stopping <- setdiff(stopping_columns, names(sims))
    if (length(missing_stopping) > 0L) {
      stop(
        "'data' is missing required stopping column(s): ",
        paste(missing_stopping, collapse = ", "),
        call. = FALSE
      )
    }
    sims$stopping_reason <- calendar_stopping_reason(
      sims$stop_futility,
      sims$stop_expected_success
    )
  }

  group_columns <- unique(c(normalized$group_vars, "scenario"))
  group_columns <- group_columns[group_columns %in% names(sims)]
  counts <- sims |>
    dplyr::group_by(!!!rlang::syms(group_columns)) |>
    dplyr::summarise(
      .calendar_n_requested = simulation_summary_denominator(
        .data$.goldilocks_requested
      ),
      .calendar_n_analyzed = dplyr::n(),
      .groups = "drop"
    )
  sims <- dplyr::left_join(sims, counts, by = group_columns)

  overall <- sims
  overall$stopping_reason <- "overall"
  duration_data <- dplyr::bind_rows(sims, overall)
  duration <- duration_data |>
    dplyr::group_by(
      !!!rlang::syms(c(group_columns, "stopping_reason"))
    ) |>
    dplyr::summarise(
      n_requested = dplyr::first(.data$.calendar_n_requested),
      n_analyzed = dplyr::first(.data$.calendar_n_analyzed),
      n_failed = .data$n_requested - .data$n_analyzed,
      n_trials = dplyr::n(),
      percent_trials = 100 * .data$n_trials / .data$n_requested,
      mean_N = mean(.data$N_enrolled),
      se_N = calendar_se(.data$N_enrolled),
      accrual_stop_mean = mean(.data$accrual_stop_time),
      accrual_stop_se = calendar_se(.data$accrual_stop_time),
      accrual_stop_p10 = calendar_quantile(.data$accrual_stop_time, 0.10),
      accrual_stop_median = calendar_quantile(.data$accrual_stop_time, 0.50),
      accrual_stop_p90 = calendar_quantile(.data$accrual_stop_time, 0.90),
      analysis_ready_mean = mean(.data$analysis_ready_time),
      analysis_ready_se = calendar_se(.data$analysis_ready_time),
      analysis_ready_p10 = calendar_quantile(.data$analysis_ready_time, 0.10),
      analysis_ready_median = calendar_quantile(.data$analysis_ready_time, 0.50),
      analysis_ready_p90 = calendar_quantile(.data$analysis_ready_time, 0.90),
      planned_completion_mean = mean(.data$planned_completion_time),
      planned_completion_se = calendar_se(.data$planned_completion_time),
      planned_completion_p10 = calendar_quantile(
        .data$planned_completion_time,
        0.10
      ),
      planned_completion_median = calendar_quantile(
        .data$planned_completion_time,
        0.50
      ),
      planned_completion_p90 = calendar_quantile(
        .data$planned_completion_time,
        0.90
      ),
      followup_person_time_mean = mean(.data$followup_person_time),
      followup_person_time_se = calendar_se(.data$followup_person_time),
      followup_person_time_p10 = calendar_quantile(
        .data$followup_person_time,
        0.10
      ),
      followup_person_time_median = calendar_quantile(
        .data$followup_person_time,
        0.50
      ),
      followup_person_time_p90 = calendar_quantile(
        .data$followup_person_time,
        0.90
      ),
      peak_active_followup_mean = mean(.data$peak_active_followup),
      peak_active_followup_se = calendar_se(.data$peak_active_followup),
      peak_active_followup_p10 = calendar_quantile(
        .data$peak_active_followup,
        0.10
      ),
      peak_active_followup_median = calendar_quantile(
        .data$peak_active_followup,
        0.50
      ),
      peak_active_followup_p90 = calendar_quantile(
        .data$peak_active_followup,
        0.90
      ),
      .groups = "drop"
    )
  duration$stopping_reason <- factor(
    duration$stopping_reason,
    levels = c(
      "expected_success",
      "futility",
      "maximum_sample_size",
      "overall"
    )
  )
  duration <- duration[do.call(order, duration[c(
    group_columns,
    "stopping_reason"
  )]), , drop = FALSE]
  duration$stopping_reason <- as.character(duration$stopping_reason)
  rownames(duration) <- NULL

  trace_input <- normalize_calendar_trace_input(data)
  if (!trace_input$available) {
    warning(
      "Interim timing is unavailable; rerun sim_trials(return_trace = TRUE). ",
      "Trial-duration summaries are still available.",
      call. = FALSE
    )
  }
  interim <- summarise_calendar_traces(trace_input$data)

  out <- list(trial_duration = duration, interim_timing = interim)
  attr(out, "time_origin") <- "first patient enrolled (time zero)"
  attr(out, "simulation_metadata") <- normalized$metadata
  class(out) <- "goldilocks_calendar_summary"
  out
}

#' Print a calendar-time operating-characteristic summary
#'
#' @param x A `goldilocks_calendar_summary` result returned by
#'   [summarise_calendar_time()].
#' @param digits A single non-negative integer giving the number of digits after
#'   the decimal point in displayed values. The default is `1`.
#' @param ... Additional arguments; currently ignored.
#'
#' @return `x`, invisibly.
#'
#' @export
print.goldilocks_calendar_summary <- function(x, digits = 1, ...) {
  cat("Calendar-time operating characteristics\n")
  cat("Time zero: first patient enrolled\n\n")

  duration <- calendar_duration_display(x$trial_duration, digits)
  cat("Trial duration and follow-up burden\n")
  print(duration, row.names = FALSE, right = FALSE)

  if (nrow(x$interim_timing) > 0L) {
    cat("\nInterim-look timing\n")
    interim <- calendar_interim_display(x$interim_timing, digits)
    print(interim, row.names = FALSE, right = FALSE)
  } else {
    cat("\nInterim-look timing: not available\n")
  }
  invisible(x)
}

#' @noRd
trial_calendar_metrics <- function(data, end_of_study) {
  enrollment <- data$enrollment
  completion <- enrollment + data$time
  accrual_stop <- max(enrollment)
  list(
    accrual_stop_time = accrual_stop,
    analysis_ready_time = max(accrual_stop, completion),
    planned_completion_time = accrual_stop + end_of_study,
    followup_person_time = sum(data$time),
    peak_active_followup = peak_active_followup(enrollment, completion)
  )
}

#' @noRd
peak_active_followup <- function(enrollment, completion) {
  times <- sort(unique(c(enrollment, completion)))
  starts <- tabulate(match(enrollment, times), nbins = length(times))
  finishes <- tabulate(match(completion, times), nbins = length(times))
  as.integer(max(c(0L, cumsum(starts - finishes))))
}

#' @noRd
active_followup_at <- function(data, calendar_time) {
  sum(
    data$enrollment <= calendar_time &
      data$enrollment + data$time > calendar_time
  )
}

#' @noRd
calendar_stopping_reason <- function(stop_futility, stop_expected_success) {
  ifelse(
    stop_expected_success != 0,
    "expected_success",
    ifelse(stop_futility != 0, "futility", "maximum_sample_size")
  )
}

#' @noRd
calendar_se <- function(x) {
  stats::sd(x) / sqrt(length(x))
}

#' @noRd
calendar_quantile <- function(x, probability) {
  unname(stats::quantile(x, probability, names = FALSE))
}

#' @noRd
normalize_calendar_trace_input <- function(data) {
  if (is.data.frame(data)) {
    return(list(data = empty_calendar_trace(), available = FALSE))
  }

  if (is_complete_sim_trials_result(data)) {
    return(normalize_calendar_trace_component(data, "1"))
  }

  scenario_names <- names(data)
  if (is.null(scenario_names)) {
    scenario_names <- as.character(seq_along(data))
  } else {
    missing_names <- is.na(scenario_names) | !nzchar(scenario_names)
    scenario_names[missing_names] <- as.character(which(missing_names))
  }
  components <- Map(normalize_calendar_trace_component, data, scenario_names)
  list(
    data = dplyr::bind_rows(lapply(components, `[[`, "data")),
    available = all(vapply(components, `[[`, logical(1), "available"))
  )
}

#' @noRd
normalize_calendar_trace_component <- function(x, scenario) {
  if (!is_complete_sim_trials_result(x) || is.null(x$traces)) {
    return(list(data = empty_calendar_trace(), available = FALSE))
  }
  metadata <- simulation_result_metadata(x, scenario)
  trace <- x$traces
  trace$scenario <- rep.int(scenario, nrow(trace))
  trace$.calendar_n_requested <- rep.int(metadata$n_requested, nrow(trace))
  trace$.calendar_n_analyzed <- rep.int(metadata$n_analyzed, nrow(trace))
  list(data = trace, available = TRUE)
}

#' @noRd
empty_calendar_trace <- function() {
  data.frame(
    scenario = character(),
    trial = integer(),
    look = integer(),
    planned_N = integer(),
    calendar_time = numeric(),
    active_followup = integer(),
    .calendar_n_requested = integer(),
    .calendar_n_analyzed = integer(),
    stringsAsFactors = FALSE
  )
}

#' @noRd
summarise_calendar_traces <- function(trace) {
  columns <- c(
    "scenario",
    "look",
    "planned_N",
    "n_requested",
    "n_analyzed",
    "n_failed",
    "n_reached",
    "percent_reached",
    "calendar_time_mean",
    "calendar_time_se",
    "calendar_time_p10",
    "calendar_time_median",
    "calendar_time_p90",
    "active_followup_mean",
    "active_followup_se",
    "active_followup_p10",
    "active_followup_median",
    "active_followup_p90"
  )
  if (nrow(trace) == 0L) {
    out <- data.frame(
      scenario = character(),
      look = integer(),
      planned_N = integer(),
      n_requested = integer(),
      n_analyzed = integer(),
      n_failed = integer(),
      n_reached = integer(),
      percent_reached = numeric(),
      calendar_time_mean = numeric(),
      calendar_time_se = numeric(),
      calendar_time_p10 = numeric(),
      calendar_time_median = numeric(),
      calendar_time_p90 = numeric(),
      active_followup_mean = numeric(),
      active_followup_se = numeric(),
      active_followup_p10 = numeric(),
      active_followup_median = numeric(),
      active_followup_p90 = numeric(),
      stringsAsFactors = FALSE
    )
    return(out[columns])
  }
  required <- c("calendar_time", "active_followup", "look", "planned_N")
  missing <- setdiff(required, names(trace))
  if (length(missing) > 0L) {
    stop(
      "Retained traces are missing required calendar-time column(s): ",
      paste(missing, collapse = ", "),
      ". Rerun the simulations with this version of goldilocks.",
      call. = FALSE
    )
  }

  out <- trace |>
    dplyr::group_by(.data$scenario, .data$look, .data$planned_N) |>
    dplyr::summarise(
      n_requested = dplyr::first(.data$.calendar_n_requested),
      n_analyzed = dplyr::first(.data$.calendar_n_analyzed),
      n_failed = .data$n_requested - .data$n_analyzed,
      n_reached = if ("trial" %in% names(trace)) {
        dplyr::n_distinct(.data$trial)
      } else {
        dplyr::n()
      },
      percent_reached = 100 * .data$n_reached / .data$n_requested,
      calendar_time_mean = mean(.data$calendar_time),
      calendar_time_se = calendar_se(.data$calendar_time),
      calendar_time_p10 = calendar_quantile(.data$calendar_time, 0.10),
      calendar_time_median = calendar_quantile(.data$calendar_time, 0.50),
      calendar_time_p90 = calendar_quantile(.data$calendar_time, 0.90),
      active_followup_mean = mean(.data$active_followup),
      active_followup_se = calendar_se(.data$active_followup),
      active_followup_p10 = calendar_quantile(.data$active_followup, 0.10),
      active_followup_median = calendar_quantile(.data$active_followup, 0.50),
      active_followup_p90 = calendar_quantile(.data$active_followup, 0.90),
      .groups = "drop"
    )
  out[columns]
}

#' @noRd
calendar_duration_display <- function(data, digits) {
  group_columns <- setdiff(
    names(data)[seq_len(match("stopping_reason", names(data)))],
    "stopping_reason"
  )
  out <- data[group_columns]
  out$stopping_reason <- gsub("_", " ", data$stopping_reason, fixed = TRUE)
  out$trials <- sprintf(
    paste0("%d (%.", digits, "f%%)"),
    data$n_trials,
    data$percent_trials
  )
  out$mean_enrolled <- round(data$mean_N, digits)
  out$accrual_stopped <- calendar_interval(
    data$accrual_stop_median,
    data$accrual_stop_p10,
    data$accrual_stop_p90,
    digits
  )
  out$analysis_ready <- calendar_interval(
    data$analysis_ready_median,
    data$analysis_ready_p10,
    data$analysis_ready_p90,
    digits
  )
  out$mean_person_time <- round(data$followup_person_time_mean, digits)
  out$mean_peak_followup <- round(data$peak_active_followup_mean, digits)
  names(out)[names(out) == "stopping_reason"] <- "stopping reason"
  names(out)[names(out) == "mean_enrolled"] <- "mean enrolled"
  names(out)[names(out) == "accrual_stopped"] <- "accrual stopped, median [P10-P90]"
  names(out)[names(out) == "analysis_ready"] <- "analysis ready, median [P10-P90]"
  names(out)[names(out) == "mean_person_time"] <- "mean person-time"
  names(out)[names(out) == "mean_peak_followup"] <- "mean peak follow-up"
  out
}

#' @noRd
calendar_interim_display <- function(data, digits) {
  out <- data[c("scenario", "look", "planned_N")]
  out$reached <- sprintf(
    paste0("%d (%.", digits, "f%%)"),
    data$n_reached,
    data$percent_reached
  )
  out$calendar_time <- calendar_interval(
    data$calendar_time_median,
    data$calendar_time_p10,
    data$calendar_time_p90,
    digits
  )
  out$active_followup <- calendar_interval(
    data$active_followup_median,
    data$active_followup_p10,
    data$active_followup_p90,
    digits
  )
  names(out)[names(out) == "planned_N"] <- "planned N"
  names(out)[names(out) == "calendar_time"] <- "calendar time, median [P10-P90]"
  names(out)[names(out) == "active_followup"] <- "active follow-up, median [P10-P90]"
  out
}

#' @noRd
calendar_interval <- function(median, lower, upper, digits) {
  format_string <- paste0("%.", digits, "f [%.", digits, "f-%.", digits, "f]")
  sprintf(format_string, median, lower, upper)
}
