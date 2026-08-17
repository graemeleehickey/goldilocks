#' @title Summarize simulations to get operating characteristics
#'
#' @param data A complete result returned by [sim_trials()], a simulation
#'   `data.frame`, or a list of either form. Named list elements identify
#'   scenarios. Existing `scenario` columns and grouping variables are
#'   preserved.
#'
#' @return Data frame reporting the operating characteristics, including the
#'   power (which will be equal to the type I error in the null case); the
#'   proportion of trials that stopped for early expected success, futility, or
#'   went to the maximum sample size. The average stopping sample size (and
#'   standard deviation) are also recorded. The proportion of trials that
#'   stopped early for expected success, yet went to ultimately fail are also
#'   reported. When complete [sim_trials()] results are supplied, the output
#'   also reports the requested, analyzed, and failed trial counts plus the
#'   effective backend and seed. Full call, RNG, parallel, timing, failure, and
#'   design metadata are retained in the `simulation_metadata` attribute.
#'
#' @importFrom dplyr bind_rows group_by summarise
#' @importFrom rlang .data
#' @export
summarise_sims <- function(data) {
  normalized <- normalize_summarise_sims_input(data)
  data <- normalized$data
  required <- c(
    "stop_futility",
    "post_prob_ha",
    "prob_threshold",
    "stop_expected_success",
    "N_enrolled"
  )
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop(
      "'data' is missing required simulation column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  group_columns <- unique(c(normalized$group_vars, "scenario"))
  group_columns <- group_columns[group_columns %in% names(data)]
  out <- data |>
    group_by(!!!rlang::syms(group_columns)) |>
    summarise(
      "power" = mean(
        !.data$stop_futility & .data$post_prob_ha > .data$prob_threshold
      ) *
        dplyr::n() /
        simulation_summary_denominator(
          .data$.goldilocks_requested
        ),
      "stop_success" = sum(.data$stop_expected_success) /
        simulation_summary_denominator(.data$.goldilocks_requested),
      "stop_futility" = sum(.data$stop_futility) /
        simulation_summary_denominator(.data$.goldilocks_requested),
      "stop_max_N" = (dplyr::n() -
        sum(.data$stop_expected_success) -
        sum(.data$stop_futility)) /
        simulation_summary_denominator(.data$.goldilocks_requested),
      "mean_N" = mean(.data$N_enrolled),
      "sd_N" = sd(.data$N_enrolled),
      "stop_and_fail" = sum(
        .data$stop_expected_success & .data$post_prob_ha <= .data$prob_threshold
      ) /
        simulation_summary_denominator(.data$.goldilocks_requested),
      ".goldilocks_n_requested" = simulation_summary_denominator(
        .data$.goldilocks_requested
      ),
      ".goldilocks_n_analyzed" = dplyr::n(),
      ".goldilocks_n_failed" = .data$.goldilocks_n_requested -
        .data$.goldilocks_n_analyzed,
      ".goldilocks_backend" = simulation_summary_metadata_value(
        .data$.goldilocks_backend,
        NA_character_
      ),
      ".goldilocks_seed" = simulation_summary_metadata_value(
        .data$.goldilocks_seed,
        NA_real_
      ),
      .groups = "drop"
    )

  internal_columns <- c(
    ".goldilocks_n_requested",
    ".goldilocks_n_analyzed",
    ".goldilocks_n_failed",
    ".goldilocks_backend",
    ".goldilocks_seed"
  )
  metric_columns <- c(
    "power",
    "stop_success",
    "stop_futility",
    "stop_max_N",
    "mean_N",
    "sd_N",
    "stop_and_fail"
  )
  if (normalized$has_sim_trials_result) {
    names(out)[match(internal_columns, names(out))] <- c(
      "n_requested",
      "n_analyzed",
      "n_failed",
      "backend",
      "seed"
    )
    out <- out[c(
      group_columns,
      "backend",
      "seed",
      "n_requested",
      "n_analyzed",
      "n_failed",
      metric_columns
    )]
    attr(out, "simulation_metadata") <- normalized$metadata

    if (!is.null(normalized$single_result)) {
      for (attribute_name in c(
        "enrollment_design",
        "decision_design",
        "rng_metadata",
        "parallel_metadata",
        "timing",
        "failures"
      )) {
        attribute_value <- attr(
          normalized$single_result,
          attribute_name,
          exact = TRUE
        )
        if (!is.null(attribute_value)) {
          attr(out, attribute_name) <- attribute_value
        }
      }
    }
  } else {
    out[internal_columns] <- NULL
    out <- out[c(group_columns, metric_columns)]
  }

  out
}

#' @title Normalize simulation-summary inputs
#'
#' @description Converts data frames, complete `sim_trials()` results, and
#'   scenario lists into one simulation data frame while retaining result
#'   metadata separately.
#'
#' @noRd
normalize_summarise_sims_input <- function(data) {
  if (is.data.frame(data)) {
    group_vars <- dplyr::group_vars(data)
    if (!"scenario" %in% names(data)) {
      data$scenario <- 1
    }
    data <- add_simulation_summary_markers(data)
    return(list(
      data = data,
      group_vars = group_vars,
      metadata = NULL,
      has_sim_trials_result = FALSE,
      single_result = NULL
    ))
  }

  if (is_complete_sim_trials_result(data)) {
    component <- normalize_simulation_summary_component(data, "1")
    return(list(
      data = component$data,
      group_vars = component$group_vars,
      metadata = stats::setNames(list(component$metadata), "1"),
      has_sim_trials_result = TRUE,
      single_result = data
    ))
  }

  if (!is.list(data) || length(data) == 0L) {
    stop(
      "'data' must be a simulation data frame, a complete sim_trials() ",
      "result, or a non-empty list of either form",
      call. = FALSE
    )
  }

  valid <- vapply(
    data,
    function(x) is.data.frame(x) || is_complete_sim_trials_result(x),
    logical(1)
  )
  if (!all(valid)) {
    stop(
      "'data' contains invalid list element(s) at position(s): ",
      paste(which(!valid), collapse = ", "),
      ". Each element must be a simulation data frame or a complete ",
      "sim_trials() result.",
      call. = FALSE
    )
  }

  scenario_names <- names(data)
  if (is.null(scenario_names)) {
    scenario_names <- as.character(seq_along(data))
  } else {
    missing_names <- is.na(scenario_names) | !nzchar(scenario_names)
    scenario_names[missing_names] <- as.character(which(missing_names))
  }

  components <- Map(
    normalize_simulation_summary_component,
    data,
    scenario_names
  )
  metadata_components <- vapply(
    components,
    function(x) !is.null(x$metadata),
    logical(1)
  )
  metadata <- lapply(components[metadata_components], `[[`, "metadata")
  names(metadata) <- scenario_names[metadata_components]

  list(
    data = bind_rows(lapply(components, `[[`, "data")),
    group_vars = unique(unlist(lapply(components, `[[`, "group_vars"))),
    metadata = metadata,
    has_sim_trials_result = any(metadata_components),
    single_result = NULL
  )
}

#' @title Normalize one simulation-summary component
#'
#' @noRd
normalize_simulation_summary_component <- function(x, scenario) {
  if (is_complete_sim_trials_result(x)) {
    metadata <- simulation_result_metadata(x, scenario)
    data <- x$sims
    data$.goldilocks_requested <- metadata$n_requested
    data$.goldilocks_backend <- metadata$backend
    data$.goldilocks_seed <- metadata$seed
  } else {
    metadata <- NULL
    data <- add_simulation_summary_markers(x)
  }
  group_vars <- setdiff(dplyr::group_vars(data), "scenario")
  data$scenario <- scenario

  list(data = data, group_vars = group_vars, metadata = metadata)
}

#' @title Add neutral simulation-summary metadata markers
#'
#' @noRd
add_simulation_summary_markers <- function(data) {
  data$.goldilocks_requested <- NA_integer_
  data$.goldilocks_backend <- NA_character_
  data$.goldilocks_seed <- NA_real_
  data
}

#' @title Identify complete sim_trials results
#'
#' @noRd
is_complete_sim_trials_result <- function(x) {
  is.list(x) &&
    is.data.frame(x$sims) &&
    is.call(x$call)
}

#' @title Extract metadata from a complete simulation result
#'
#' @noRd
simulation_result_metadata <- function(x, scenario) {
  parallel_metadata <- attr(x, "parallel_metadata", exact = TRUE)
  rng_metadata <- attr(x, "rng_metadata", exact = TRUE)
  failures <- if (!is.null(x$failures)) {
    x$failures
  } else {
    attr(x, "failures", exact = TRUE)
  }
  known_failures <- simulation_failure_count(failures)
  call_trials <- x$call$N_trials
  if (!is.numeric(call_trials) || length(call_trials) != 1L) {
    call_trials <- NULL
  }
  requested_candidates <- c(
    nrow(x$sims) + known_failures,
    parallel_metadata$tasks,
    call_trials
  )
  requested_candidates <- requested_candidates[
    is.finite(requested_candidates) & requested_candidates >= nrow(x$sims)
  ]
  n_requested <- if (length(requested_candidates) == 0L) {
    nrow(x$sims)
  } else {
    max(requested_candidates)
  }
  backend <- parallel_metadata$backend
  if (is.null(backend)) {
    backend <- rng_metadata$backend
  }
  if (is.null(backend)) {
    backend <- NA_character_
  }
  seed <- rng_metadata$stream_seed
  if (is.null(seed)) {
    call_seed <- x$call$seed
    seed <- if (is.numeric(call_seed) && length(call_seed) == 1L) {
      call_seed
    } else {
      NA_real_
    }
  }

  list(
    scenario = scenario,
    backend = as.character(backend)[1],
    seed = as.numeric(seed)[1],
    n_requested = as.integer(n_requested),
    n_analyzed = as.integer(nrow(x$sims)),
    n_failed = as.integer(n_requested - nrow(x$sims)),
    call = x$call,
    enrollment_design = attr(x, "enrollment_design", exact = TRUE),
    decision_design = attr(x, "decision_design", exact = TRUE),
    rng_metadata = rng_metadata,
    parallel_metadata = parallel_metadata,
    timing = if (!is.null(x$timing)) {
      x$timing
    } else {
      attr(x, "timing", exact = TRUE)
    },
    failures = failures
  )
}

#' @title Count recorded simulation failures
#'
#' @noRd
simulation_failure_count <- function(failures) {
  if (is.null(failures)) {
    return(0L)
  }
  if (is.numeric(failures) && length(failures) == 1L) {
    return(as.integer(failures))
  }
  if (is.data.frame(failures)) {
    return(as.integer(nrow(failures)))
  }
  as.integer(length(failures))
}

#' @title Determine a simulation-summary denominator
#'
#' @noRd
simulation_summary_denominator <- function(requested) {
  requested <- requested[!is.na(requested)]
  if (length(requested) == 0L) {
    return(dplyr::n())
  }
  max(requested)
}

#' @title Collapse scalar simulation metadata within a scenario
#'
#' @noRd
simulation_summary_metadata_value <- function(x, missing) {
  x <- unique(x[!is.na(x)])
  if (length(x) == 0L) missing else x[1]
}
