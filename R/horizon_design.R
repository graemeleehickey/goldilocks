# Validate a prespecified schedule against the looks and horizons in saved paths.
validate_endpoint_schedule <- function(endpoint_schedule, args) {
  if (
    !is.data.frame(endpoint_schedule) ||
      !identical(names(endpoint_schedule), c("N", "horizon")) ||
      !is.numeric(endpoint_schedule$N) ||
      !is.numeric(endpoint_schedule$horizon) ||
      anyNA(endpoint_schedule) ||
      !identical(
        as.numeric(endpoint_schedule$N),
        as.numeric(args$interim_look)
      ) ||
      any(!endpoint_schedule$horizon %in% args$endpoint_horizons) ||
      any(diff(endpoint_schedule$horizon) < 0)
  ) {
    stop(
      "'endpoint_schedule' must have columns N and horizon, with every planned interim N in order and nondecreasing horizons available in the saved paths"
    )
  }
  endpoint_schedule
}

validate_horizon_arguments <- function(args) {
  if ("endpoint_schedule" %in% names(args)) {
    schedule <- args$endpoint_schedule
    if (
      !is.data.frame(schedule) ||
        !identical(names(schedule), c("N", "horizon")) ||
        !is.numeric(schedule$N) ||
        !is.numeric(schedule$horizon) ||
        any(!is.finite(as.matrix(schedule))) ||
        any(schedule$horizon <= 0) ||
        any(diff(schedule$horizon) < 0)
    ) {
      stop(
        "'endpoint_schedule' must have numeric columns N and horizon, finite values, and positive, nondecreasing horizons"
      )
    }
    args["interim_look"] <- list(if (nrow(schedule)) schedule$N else NULL)
    args$endpoint_horizons <- sort(unique(c(schedule$horizon, args$final_horizon)))
    args$followup_horizon <- args$final_horizon
  }
  # Older saved paths retain their original grid and collection duration.
  horizons <- args$endpoint_horizons
  if (
    !is.numeric(horizons) ||
      !is.null(dim(horizons)) ||
      !length(horizons) ||
      any(!is.finite(horizons)) ||
      any(horizons <= 0) ||
      any(diff(horizons) <= 0)
  ) {
    stop(
      "'endpoint_horizons' must be positive, finite, and strictly increasing"
    )
  }
  for (name in c("final_horizon", "followup_horizon", "dropout_horizon")) {
    validate_endpoint_time(args[[name]], NULL, name)
  }
  if (
    max(horizons) != args$final_horizon ||
      args$followup_horizon < args$final_horizon
  ) {
    stop(
      "'final_horizon' must be at least every scheduled horizon, and saved collection must cover it"
    )
  }
  args$backend <- match.arg(
    args$backend,
    c("auto", "fork", "psock", "sequential")
  )
  args$empty_interval <- match.arg(
    args$empty_interval,
    c("prior", "propagate", "error")
  )
  args$end_of_study <- args$followup_horizon
  args$binary_imputation <- "event-time"
  args$rmst_tau <- args$final_horizon
  args <- validate_path_arguments(args)

  settings <- args$horizon_settings
  if (is.null(settings)) {
    settings <- data.frame(horizon = horizons)
  }
  if (
    !is.data.frame(settings) ||
      !"horizon" %in% names(settings) ||
      !is.numeric(settings$horizon) ||
      anyNA(settings$horizon) ||
      anyDuplicated(names(settings)) ||
      anyDuplicated(settings$horizon) ||
      !setequal(settings$horizon, horizons) ||
      any(!names(settings) %in% c("horizon", "h0", "rmst_tau", "prior_bin"))
  ) {
    stop(
      "'horizon_settings' must have one row per endpoint horizon and columns horizon, optionally h0, rmst_tau, and prior_bin"
    )
  }
  settings <- settings[match(horizons, settings$horizon), , drop = FALSE]
  if (!"h0" %in% names(settings)) {
    settings$h0 <- args$h0
  }
  if (!"rmst_tau" %in% names(settings)) {
    settings$rmst_tau <- horizons
  }
  if (!"prior_bin" %in% names(settings)) {
    settings$prior_bin <- rep(list(args$prior_bin), length(horizons))
  }
  if (!is.list(settings$prior_bin)) {
    stop("'horizon_settings$prior_bin' must be a list column")
  }
  for (i in seq_along(horizons)) {
    validate_h0(settings$h0[i], args$method, is.null(args$hazard_control))
    validate_rmst_args(
      settings$rmst_tau[i],
      horizons[i],
      if (args$method == "rmst") settings$h0[i] else 0
    )
    if (args$method == "bayes-bin") {
      validate_bayes_binomial_args(
        settings$prior_bin[[i]],
        args$bin_method,
        args$N_mcmc
      )
    }
  }
  args$horizon_settings <- settings
  args
}

# Derive each endpoint analysis from the master grid without moving late priors
# into earlier intervals. A cutpoint equal to the endpoint starts no new interval.
horizon_analysis_config <- function(args, horizon) {
  config <- args
  config$single_arm <- is.null(args$hazard_control)
  config$binary_imputation <- "event-time"
  config$end_of_study <- horizon
  cutpoints <- args$cutpoints[args$cutpoints < horizon]
  config["cutpoints"] <- list(if (length(cutpoints)) cutpoints else NULL)
  intervals <- seq_len(length(config$cutpoints) + 1L)
  for (name in c("prior_surv", "prior_surv_final")) {
    config[[name]] <- args[[name]][, intervals, , drop = FALSE]
  }
  row <- match(horizon, args$horizon_settings$horizon)
  config$h0 <- args$horizon_settings$h0[row]
  config$rmst_tau <- args$horizon_settings$rmst_tau[row]
  config$prior_bin <- args$horizon_settings$prior_bin[[row]]
  if (!uses_posterior_draws(config$method, config$bin_method)) {
    config$N_mcmc <- 1L
  }
  config
}

# The generated histories share one event and dropout process. Re-censoring
# changes endpoint status, but does not redraw either process.
censor_horizon_data <- function(data, horizon) {
  data$event <- as.numeric(data$event == 1 & data$time <= horizon)
  data$loss_to_fu <- data$loss_to_fu & data$time < horizon
  data$time <- pmin(data$time, horizon)
  data
}

horizon_interim_evidence <- function(data, config, look, n) {
  cut <- censor_horizon_data(data, config$end_of_study)
  call_path_function(
    evaluate_interim_evidence,
    config,
    data_interim = prepare_simulated_interim(cut, n, config$end_of_study),
    look = look,
    planned_N = n,
    calendar_time = data$enrollment[n],
    active_followup = active_followup_at(cut, data$enrollment[n]),
    check_futility = config$end_of_study == config$final_horizon,
    isolate_errors = TRUE
  )
}

# Keep the current-cohort evidence at its own data cut. Only the maximum-cohort
# scores come from the final-horizon fit, even when the earlier endpoint stops.
combine_horizon_evidence <- function(current, maximum, horizon, final_horizon) {
  if (is.null(current)) {
    return(NULL)
  }
  current$scores$maximum <- if (is.null(maximum)) {
    structure(
      numeric(),
      failure = list(message = "Final-horizon prediction failed")
    )
  } else {
    maximum$scores$maximum
  }
  current$context$success_horizon <- horizon
  current$context$futility_horizon <- final_horizon
  current$context$success_data_horizon <- horizon
  current$context$futility_data_horizon <- final_horizon
  for (name in c("events_treatment", "events_control", "N_pending")) {
    current$context[[paste0("futility_", name)]] <- if (is.null(maximum)) {
      NA_integer_
    } else {
      maximum$context[[name]]
    }
  }
  for (name in c("warnings", "empty_interval_fallbacks")) {
    current$diagnostics[[name]] <- unique(c(
      current$diagnostics[[name]],
      maximum$diagnostics[[name]]
    ))
  }
  current
}

normalize_horizon_rules <- function(args, Fn, Sn, Qn, prob_ha) {
  validate_single_probability(prob_ha, "prob_ha")
  n <- length(args$interim_look)
  rules <- list(
    Fn = normalize_interim_threshold(Fn, n, "Fn", null_disables = TRUE),
    Sn = normalize_interim_threshold(Sn, n, "Sn"),
    Qn = normalize_interim_threshold(Qn, n, "Qn"),
    prob_ha = prob_ha
  )
  validate_success_threshold_order(rules$Sn, rules$Qn)
  rules
}

horizon_interim_decision <- function(evidence, rules, i) {
  if (is.null(evidence)) {
    stop("Required success prediction failed")
  }
  require_path_scores(evidence$scores$current, "Success prediction")
  now <- mean(evidence$scores$current > rules$prob_ha)
  maximum <- NA_real_
  if (now <= rules$Sn[i] && rules$Fn[i] > 0) {
    require_path_scores(evidence$scores$maximum, "Futility prediction")
    maximum <- mean(evidence$scores$maximum > rules$prob_ha)
  }
  select_interim_decision(
    now,
    maximum,
    rules$Sn[i],
    rules$Qn[i],
    rules$Fn[i],
    rules$Fn[i] > 0
  )$decision
}
