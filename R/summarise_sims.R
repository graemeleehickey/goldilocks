#' @title Summarize simulations to get operating characteristics
#'
#' @param data A complete result returned by [sim_trials()], a simulation
#'   `data.frame`, or a list of either form. Named list elements identify
#'   scenarios. Existing `scenario` columns and grouping variables are
#'   preserved.
#' @param max_mcse Optional named numeric vector giving the largest acceptable
#'   Monte Carlo standard error for selected estimands. Supported names are
#'   `power`, `stop_success`, `stop_futility`, `stop_max_N`, `mean_N`,
#'   `stop_and_fail`, and `failure_rate`. A warning identifies every scenario
#'   and estimand whose Monte Carlo standard error exceeds its target.
#'
#' @return Data frame reporting the operating characteristics, including the
#'   power (which will be equal to the type I error in the null case); the
#'   proportion of trials that stopped for early expected success, futility, or
#'   went to the maximum sample size. The average stopping sample size (and
#'   standard deviation) are also recorded. The proportion of trials that
#'   stopped early for expected success, yet went to ultimately fail are also
#'   reported. Each probability and mean is accompanied by its Monte Carlo
#'   standard error and 95% Monte Carlo confidence limits, with columns ending
#'   in `_mcse`, `_mc_lower`, and `_mc_upper`. Probability intervals use the
#'   Wilson method; the mean sample-size interval uses a t distribution.
#'
#'   These intervals describe uncertainty from using a finite number of
#'   simulated trials under fixed design and data-generating assumptions. They
#'   are **not clinical confidence intervals**, treatment-effect intervals, or
#'   measures of model uncertainty.
#'
#'   The output always reports `n_used`, the number of successfully analyzed
#'   simulations used by the operating-characteristic estimands. When complete
#'   [sim_trials()] results are supplied, it also reports requested, analyzed,
#'   and failed counts, the failure rate, effective backend, and seed. For a raw
#'   simulation data frame, requested and failed counts are unknown and are
#'   reported as `NA`. Full call, RNG, parallel, timing, failure, and design
#'   metadata, including the evaluated simulation arguments, are retained in
#'   the `simulation_metadata` attribute.
#'
#' @importFrom dplyr bind_rows group_by summarise
#' @importFrom rlang .data
#' @export
summarise_sims <- function(data, max_mcse = NULL) {
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

  probability_metrics <- c(
    "power",
    "stop_success",
    "stop_futility",
    "stop_max_N",
    "stop_and_fail"
  )
  precision_metrics <- c(probability_metrics, "mean_N", "failure_rate")
  validate_max_mcse(max_mcse, precision_metrics)

  data$.goldilocks_power <- with(
    data,
    !stop_futility & post_prob_ha > prob_threshold
  )
  data$.goldilocks_stop_success <- data$stop_expected_success != 0
  data$.goldilocks_stop_futility <- data$stop_futility != 0
  data$.goldilocks_stop_max_N <- with(
    data,
    !.goldilocks_stop_success & !.goldilocks_stop_futility
  )
  data$.goldilocks_stop_and_fail <- with(
    data,
    .goldilocks_stop_success & post_prob_ha <= prob_threshold
  )

  group_columns <- unique(c(normalized$group_vars, "scenario"))
  group_columns <- group_columns[group_columns %in% names(data)]
  out <- data |>
    group_by(!!!rlang::syms(group_columns)) |>
    summarise(
      "power" = mean(.data$.goldilocks_power),
      "stop_success" = mean(.data$.goldilocks_stop_success),
      "stop_futility" = mean(.data$.goldilocks_stop_futility),
      "stop_max_N" = mean(.data$.goldilocks_stop_max_N),
      "mean_N" = mean(.data$N_enrolled),
      "sd_N" = sd(.data$N_enrolled),
      "stop_and_fail" = mean(.data$.goldilocks_stop_and_fail),
      ".goldilocks_n_requested" = simulation_summary_requested(
        .data$.goldilocks_requested
      ),
      ".goldilocks_n_used" = dplyr::n(),
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

  out$.goldilocks_n_analyzed <- out$.goldilocks_n_used
  out$.goldilocks_n_failed <- ifelse(
    is.na(out$.goldilocks_n_requested),
    NA_integer_,
    out$.goldilocks_n_requested - out$.goldilocks_n_analyzed
  )
  out$failure_rate <- ifelse(
    is.na(out$.goldilocks_n_requested),
    NA_real_,
    out$.goldilocks_n_failed / out$.goldilocks_n_requested
  )

  for (metric in probability_metrics) {
    out <- add_binomial_mc_columns(
      out,
      metric,
      out$.goldilocks_n_used
    )
  }
  out <- add_mean_mc_columns(
    out,
    metric = "mean_N",
    standard_deviation = out$sd_N,
    denominator = out$.goldilocks_n_used
  )
  out <- add_binomial_mc_columns(
    out,
    metric = "failure_rate",
    denominator = out$.goldilocks_n_requested
  )

  internal_columns <- c(
    ".goldilocks_n_requested",
    ".goldilocks_n_analyzed",
    ".goldilocks_n_failed",
    ".goldilocks_n_used",
    ".goldilocks_backend",
    ".goldilocks_seed"
  )
  names(out)[match(internal_columns[1:4], names(out))] <- c(
    "n_requested",
    "n_analyzed",
    "n_failed",
    "n_used"
  )
  denominator_columns <- c(
    "n_requested",
    "n_analyzed",
    "n_failed",
    "n_used",
    "failure_rate",
    "failure_rate_mcse",
    "failure_rate_mc_lower",
    "failure_rate_mc_upper"
  )
  stopping_probability_metrics <- setdiff(
    probability_metrics,
    "stop_and_fail"
  )
  mean_columns <- c(
    "mean_N",
    "mean_N_mcse",
    "mean_N_mc_lower",
    "mean_N_mc_upper",
    "sd_N"
  )
  metric_columns <- c(
    unlist(lapply(
      stopping_probability_metrics,
      simulation_mc_column_names
    )),
    mean_columns,
    simulation_mc_column_names("stop_and_fail")
  )
  if (normalized$has_sim_trials_result) {
    names(out)[match(internal_columns[5:6], names(out))] <- c(
      "backend",
      "seed"
    )
    out <- out[c(
      group_columns,
      "backend",
      "seed",
      denominator_columns,
      metric_columns
    )]
    attr(out, "simulation_metadata") <- normalized$metadata

    if (!is.null(normalized$single_result)) {
      for (attribute_name in c(
        "arguments",
        "enrollment_design",
        "decision_design",
        "prior_design",
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
    out[internal_columns[5:6]] <- NULL
    out <- out[c(group_columns, denominator_columns, metric_columns)]
  }

  attr(out, "mc_conf_level") <- 0.95
  warn_simulation_precision(out, max_mcse, group_columns)

  out
}

#' @title Names used for one Monte Carlo interval
#'
#' @noRd
simulation_mc_column_names <- function(metric) {
  c(
    metric,
    paste0(metric, "_mcse"),
    paste0(metric, "_mc_lower"),
    paste0(metric, "_mc_upper")
  )
}

#' @title Add binomial Monte Carlo uncertainty columns
#'
#' @noRd
add_binomial_mc_columns <- function(data, metric, denominator) {
  uncertainty <- binomial_mc_uncertainty(data[[metric]], denominator)
  data[[paste0(metric, "_mcse")]] <- uncertainty$mcse
  data[[paste0(metric, "_mc_lower")]] <- uncertainty$lower
  data[[paste0(metric, "_mc_upper")]] <- uncertainty$upper
  data
}

#' @title Calculate binomial Monte Carlo uncertainty
#'
#' @description Uses the Wilson interval so zero and one observed outcomes do
#'   not produce misleading zero-width intervals.
#'
#' @noRd
binomial_mc_uncertainty <- function(estimate, denominator) {
  estimate <- as.numeric(estimate)
  denominator <- as.numeric(denominator)
  out <- data.frame(
    mcse = rep(NA_real_, length(estimate)),
    lower = rep(NA_real_, length(estimate)),
    upper = rep(NA_real_, length(estimate))
  )
  valid <- is.finite(estimate) &
    estimate >= 0 &
    estimate <= 1 &
    is.finite(denominator) &
    denominator > 0
  if (!any(valid)) {
    return(out)
  }

  probability <- estimate[valid]
  n <- denominator[valid]
  z <- stats::qnorm(0.975)
  adjustment <- 1 + z^2 / n
  center <- (probability + z^2 / (2 * n)) / adjustment
  half_width <- z / adjustment * sqrt(
    probability * (1 - probability) / n + z^2 / (4 * n^2)
  )
  out$mcse[valid] <- sqrt(probability * (1 - probability) / n)
  out$lower[valid] <- pmax(0, center - half_width)
  out$upper[valid] <- pmin(1, center + half_width)
  out
}

#' @title Add mean Monte Carlo uncertainty columns
#'
#' @noRd
add_mean_mc_columns <- function(
  data,
  metric,
  standard_deviation,
  denominator
) {
  uncertainty <- mean_mc_uncertainty(
    data[[metric]],
    standard_deviation,
    denominator
  )
  data[[paste0(metric, "_mcse")]] <- uncertainty$mcse
  data[[paste0(metric, "_mc_lower")]] <- uncertainty$lower
  data[[paste0(metric, "_mc_upper")]] <- uncertainty$upper
  data
}

#' @title Calculate Monte Carlo uncertainty for a mean
#'
#' @noRd
mean_mc_uncertainty <- function(estimate, standard_deviation, denominator) {
  estimate <- as.numeric(estimate)
  standard_deviation <- as.numeric(standard_deviation)
  denominator <- as.numeric(denominator)
  mcse <- standard_deviation / sqrt(denominator)
  critical_value <- rep(NA_real_, length(denominator))
  valid_denominator <- is.finite(denominator) & denominator > 1
  critical_value[valid_denominator] <- stats::qt(
    0.975,
    df = denominator[valid_denominator] - 1
  )
  half_width <- critical_value * mcse
  valid <- is.finite(estimate) &
    is.finite(mcse) &
    is.finite(critical_value) &
    valid_denominator
  data.frame(
    mcse = ifelse(valid, mcse, NA_real_),
    lower = ifelse(valid, estimate - half_width, NA_real_),
    upper = ifelse(valid, estimate + half_width, NA_real_)
  )
}

#' @title Validate requested Monte Carlo precision targets
#'
#' @noRd
validate_max_mcse <- function(max_mcse, supported) {
  if (is.null(max_mcse)) {
    return(invisible(NULL))
  }
  if (
    !is.numeric(max_mcse) ||
      length(max_mcse) == 0L ||
      any(!is.finite(max_mcse)) ||
      any(max_mcse <= 0) ||
      is.null(names(max_mcse)) ||
      any(is.na(names(max_mcse))) ||
      any(!nzchar(names(max_mcse))) ||
      anyDuplicated(names(max_mcse))
  ) {
    stop(
      "'max_mcse' must be a named numeric vector of unique, positive, ",
      "finite targets",
      call. = FALSE
    )
  }
  unknown <- setdiff(names(max_mcse), supported)
  if (length(unknown) > 0L) {
    stop(
      "Unsupported 'max_mcse' estimand(s): ",
      paste(unknown, collapse = ", "),
      ". Supported estimands are: ",
      paste(supported, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' @title Warn when Monte Carlo precision targets are not met
#'
#' @noRd
warn_simulation_precision <- function(data, max_mcse, group_columns) {
  if (is.null(max_mcse)) {
    return(invisible(NULL))
  }
  scenario_labels <- apply(
    as.data.frame(data[group_columns]),
    1,
    function(values) {
      paste(paste(group_columns, values, sep = "="), collapse = ", ")
    }
  )
  problems <- character()
  for (metric in names(max_mcse)) {
    observed <- data[[paste0(metric, "_mcse")]]
    failed <- is.na(observed) | observed > max_mcse[[metric]]
    if (!any(failed)) {
      next
    }
    for (row in which(failed)) {
      detail <- if (is.na(observed[row])) {
        "MCSE is unavailable"
      } else {
        paste0(
          "MCSE ",
          formatC(observed[row], digits = 4, format = "f"),
          " exceeds ",
          formatC(max_mcse[[metric]], digits = 4, format = "f")
        )
      }
      problems <- c(
        problems,
        paste0("  - ", scenario_labels[row], ", ", metric, ": ", detail)
      )
    }
  }
  if (length(problems) > 0L) {
    warning(
      paste(
        c("Monte Carlo precision target not met:", problems),
        collapse = "\n"
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
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
    arguments = attr(x, "arguments", exact = TRUE),
    enrollment_design = attr(x, "enrollment_design", exact = TRUE),
    decision_design = attr(x, "decision_design", exact = TRUE),
    prior_design = attr(x, "prior_design", exact = TRUE),
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

#' @title Return a known requested simulation count
#'
#' @description Unlike `simulation_summary_denominator()`, this helper retains
#'   `NA` when a raw simulation data frame has discarded the requested count.
#'
#' @noRd
simulation_summary_requested <- function(requested) {
  requested <- requested[!is.na(requested)]
  if (length(requested) == 0L) {
    return(NA_integer_)
  }
  as.integer(max(requested))
}

#' @title Collapse scalar simulation metadata within a scenario
#'
#' @noRd
simulation_summary_metadata_value <- function(x, missing) {
  x <- unique(x[!is.na(x)])
  if (length(x) == 0L) missing else x[1]
}
