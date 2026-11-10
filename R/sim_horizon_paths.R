#' Simulate complete trial paths for prespecified endpoint horizons
#'
#' @description Generates one set of event and dropout histories per trial and evaluates every planned interim at each supported endpoint horizon. The saved evidence can be reused by [apply_horizon_rules()] to calibrate stopping thresholds and prespecified horizon schedules.
#' @inheritParams sim_trial_paths
#' @param endpoint_schedule Data frame with exactly two numeric columns, `N` and `horizon`, containing every planned interim enrollment size in increasing order and a positive, nondecreasing endpoint horizon at each look. All interim sizes must be less than `N_total`. The saved horizons are derived as `sort(unique(c(endpoint_schedule$horizon, final_horizon)))` and evaluated at every look, allowing subsequent comparison of schedules using the same looks and endpoint times. Use a zero-row data frame with numeric columns when there are no interims.
#' @param final_horizon Required endpoint time for the maximum-\eqn{N} final analysis, all futility predictions, and planned outcome collection for every participant from enrollment. It must be at least every scheduled horizon. The generating hazards must describe this entire period.
#' @param dropout_horizon Reference time at which `prop_loss` specifies the independent exponential dropout probability. Defaults to `final_horizon`; it remains fixed when the analysis endpoint changes.
#' @param prop_loss Dropout probability by `dropout_horizon`, with the same arm-specific conventions as [sim_comp_data()]. The corresponding dropout hazard is held fixed across all endpoint horizons.
#' @param horizon_settings Optional data frame with one row per distinct scheduled horizon plus `final_horizon`, identified by column `horizon`. Optional columns are `h0`, `rmst_tau`, and a list column `prior_bin`. Missing columns use the common `h0`, the endpoint horizon itself for RMST restriction, and the common `prior_bin`, respectively. These settings are fixed within a set of trial paths.
#' @param cutpoints Master piecewise-exponential analysis cutpoints, strictly inside `final_horizon`. Each endpoint uses only cutpoints strictly before that endpoint and the corresponding leading intervals of `prior_surv` and `prior_surv_final`. If an endpoint H lies inside a master interval, the same hazard parameter and unscaled Gamma prior apply, but its success update uses events and exposure only through H. Available observations beyond H are excluded even within that interval; H is not a new hazard-change assumption in the generating model.
#' @param prior_surv_final Gamma prior for final-stage imputation and Bayesian survival analysis of completed predictive and final data. Accepts the same master-grid forms as `prior_surv` and defaults to it. At each endpoint, only the corresponding leading intervals are used. Supply both priors explicitly to separate predictive and analysis assumptions.
#' @param h0 Common null value or margin, overridden by `horizon_settings$h0` when supplied. Bayesian and risk-difference methods use event-probability scales, Cox uses log hazard ratio, RMST uses treatment-control differences in time units, and log-rank requires zero. See [sim_trials()] for method-specific directions and allowable ranges.
#' @param imputed_final Logical; whether to impute outcomes lost before the selected endpoint in final analyses. Defaults to `FALSE`. Final-imputation methods and restrictions are the same as in [sim_trials()]: log-rank imputation is unavailable, and Farrington-Manning imputation requires no missing outcomes. Cox, RMST, and risk-difference Wald imputation require at least two imputations and positive pooled variance.
#'
#' @details Success prediction at horizon H fits its predictive hazard model using only follow-up through H that is available at the look's calendar date. It imputes completed outcomes through H among currently enrolled participants. Futility prediction uses a separate model fitted to available follow-up through `final_horizon`, and predicts success at that horizon at maximum \eqn{N}. Follow-up after H cannot inform the success model, including when a hazard interval crosses H. When H equals `final_horizon`, the two predictions share a fitted predictive model.
#'
#' Event times are generated once through `final_horizon` and re-censored for each analysis. Predictive imputation uses event times; independently imputed binary outcomes at different horizons are not supported. When the success and futility endpoints differ, the two models make separate posterior and imputation draws. The first completed data set from each model is not the same future trial observed for different durations; draw numbers cannot be matched to estimate joint predictive outcomes.
#'
#' The final analysis after an expected-success stop uses that look's endpoint, even if it fails. It never resumes the trial at a longer endpoint. Futility targets the maximum-\eqn{N} final analysis, rather than success at any remaining interim. Interim sample sizes must be less than `N_total`; maximum \eqn{N} is followed by final endpoint ascertainment.
#'
#' `prop_loss` describes an exponential dropout distribution independently of the analysis endpoint. For example, `prop_loss = 0.05, dropout_horizon = 12, final_horizon = 18` retains a 5% dropout probability by month 12 while generating follow-up through month 18. Dropout after month 12 does not make the 12-month endpoint missing.
#'
#' All other fixed design settings and analysis restrictions follow [sim_trial_paths()]. Explicit seeds preserve the caller's random-number state and reproduce results across supported computation methods and consecutive `trial_offset` batches. Changing the set of retained horizons changes the allocation of random numbers. Changing only a candidate schedule or thresholds does not. Select designs using one set of trial paths and assess the complete adaptive design using independent simulations.
#'
#' @return A `goldilocks_horizon_paths` object with `trials`, a calculation-level `failures` table including `horizon`, `schema_version`, and `call`. Each trial retains per-look evidence for each horizon and final analyses for each cohort and horizon. Separate posterior diagnostics record the event counts, exposure, and Gamma updates used at each data cut. Attributes retain evaluated arguments, horizon-specific priors, enrollment design, and random-number and computation settings. Save and restore paths with [saveRDS()] and [readRDS()]. Memory grows with trials, looks, horizons, and imputation count.
#' @seealso [sim_horizon_trials()], [apply_horizon_rules()]
#' @export
#' @examples
#' schedule <- data.frame(N = c(20, 40), horizon = c(12, 18))
#' paths <- sim_horizon_paths(
#'   hazard_treatment = c(0.04, 0.03), hazard_control = c(0.07, 0.06),
#'   generation_cutpoints = 12, cutpoints = 12,
#'   N_total = 60, endpoint_schedule = schedule, final_horizon = 18,
#'   lambda = 3, method = "bayes-bin", bin_method = "quadrature",
#'   alternative = "less", N_impute = 5, N_trials = 2, seed = 123
#' )
#' result <- apply_horizon_rules(paths, schedule, return_trace = TRUE)
#' summarise_sims(result)
sim_horizon_paths <- function(
  hazard_treatment,
  hazard_control = NULL,
  cutpoints = NULL,
  N_total,
  lambda = 0.3,
  lambda_time = NULL,
  endpoint_schedule,
  final_horizon,
  dropout_horizon = final_horizon,
  horizon_settings = NULL,
  prior_surv = c(0.1, 0.1),
  prior_bin = c(1, 1),
  bin_method = "mc",
  block = 2,
  rand_ratio = c(control = 1, treatment = 1),
  prop_loss = 0,
  alternative = "greater",
  h0 = 0,
  N_impute = 500,
  N_mcmc = 1000,
  mc_conf_level = 0.95,
  N_trials = 10,
  method = "logrank",
  imputed_final = FALSE,
  empty_interval = c("prior", "propagate", "error"),
  ncores = 1L,
  backend = c("auto", "fork", "psock", "sequential"),
  seed = NULL,
  prior_surv_final = prior_surv,
  generation_cutpoints = cutpoints,
  trial_offset = 0L
) {
  Call <- match.call()
  if (is.function(Call[[1L]])) {
    Call[[1L]] <- quote(sim_horizon_paths)
  }
  args <- validate_horizon_arguments(capture_arguments(
    sim_horizon_paths,
    environment()
  ))
  run_horizon_simulations(args, Call)
}

run_horizon_simulations <- function(args, Call, schedule = NULL, rules = NULL) {
  execution <- resolve_sim_execution(args$backend, args$ncores, args$N_trials)
  inner_mc <- args$method == "bayes-surv" ||
    (args$method == "bayes-bin" && args$bin_method == "mc")
  score_bytes <- as.double(args$N_trials) *
    length(args$interim_look) *
    args$N_impute *
    (length(args$endpoint_horizons) + 1L) *
    if (inner_mc) 12 else 8
  if (is.null(schedule) && score_bytes >= 512 * 1024^2) {
    message(sprintf(
      "Predictive scores require approximately %.1f MiB before diagnostics and object overhead; consider separate batches with 'trial_offset'.",
      score_bytes / 1024^2
    ))
  }
  caller_kind <- RNGkind()
  stream_seed <- if (is.null(args$seed)) {
    sample.int(.Machine$integer.max - 1L, 1L)
  } else {
    args$seed
  }
  restore_rng <- preserve_path_rng()
  on.exit(restore_rng(), add = TRUE)
  trial_ids <- as.integer(args$trial_offset + seq_len(args$N_trials))
  streams <- make_rng_streams(stream_seed, max(trial_ids))[trial_ids]
  configs <- lapply(args$endpoint_horizons, function(h) {
    horizon_analysis_config(args, h)
  })
  simulate <- if (execution$backend == "psock") {
    make_psock_callable("simulate_horizon_path")
  } else {
    simulate_horizon_path
  }
  worker <- function(i) {
    simulate(args, configs, streams[[i]], trial_ids[i], schedule, rules)
  }
  indices <- seq_len(args$N_trials)
  trials <- switch(
    execution$backend,
    sequential = lapply(indices, worker),
    fork = pbmcapply::pbmclapply(indices, worker, mc.cores = execution$workers),
    psock = {
      cluster <- make_sim_cluster(execution$workers)
      on.exit(stop_sim_cluster(cluster), add = TRUE)
      initialize_sim_cluster(cluster)
      run_sim_cluster(cluster, indices, worker)
    }
  )
  failures <- dplyr::bind_rows(lapply(trials, `[[`, "failures"))
  trials <- lapply(trials, function(x) {
    x$failures <- NULL
    x
  })
  out <- structure(
    list(
      trials = trials,
      failures = failures,
      schema_version = 2L,
      call = Call
    ),
    class = "goldilocks_horizon_paths"
  )
  attr(out, "arguments") <- args[names(formals(sim_horizon_paths))]
  attr(out, "score_storage_bytes") <- score_bytes
  attr(out, "enrollment_design") <- call_path_function(
    new_enrollment_design,
    args
  )
  attr(out, "prior_design") <- dplyr::bind_rows(lapply(
    configs,
    function(config) {
      priors <- dplyr::bind_rows(
        call_path_function(gamma_prior_diagnostics, config, stage = "interim"),
        call_path_function(
          gamma_prior_diagnostics,
          config,
          prior_surv = config$prior_surv_final,
          stage = "final"
        )
      )
      priors$horizon <- config$end_of_study
      priors
    }
  ))
  attr(out, "runtime_metadata") <- list(
    package_version = as.character(utils::packageVersion("goldilocks")),
    R_version = as.character(getRversion())
  )
  attr(out, "rng_metadata") <- list(
    caller_kind = caller_kind,
    stream_kind = "L'Ecuyer-CMRG",
    seed_policy = if (is.null(args$seed)) {
      "caller_derived_bank"
    } else {
      "explicit_preserve_caller"
    },
    stream_seed = stream_seed,
    stream_layout = "trial-look-horizon-stage-v1",
    backend = execution$backend,
    ncores = as.integer(execution$workers)
  )
  attr(out, "parallel_metadata") <- list(
    requested_backend = args$backend,
    backend = execution$backend,
    selection_reason = execution$reason,
    requested_ncores = as.integer(args$ncores),
    workers = as.integer(execution$workers),
    tasks = as.integer(args$N_trials)
  )
  if (nrow(failures)) {
    warning(
      nrow(failures),
      " horizon calculations failed; see 'paths$failures'. Candidate rules determine which trials are evaluable.",
      call. = FALSE
    )
  }
  out
}

simulate_horizon_path <- function(
  args,
  configs,
  stream,
  trial,
  schedule = NULL,
  rules = NULL
) {
  looks <- args$interim_look
  cohorts <- c(looks, args$N_total)
  horizons <- args$endpoint_horizons
  final_index <- match(args$final_horizon, horizons)
  failures <- empty_path_failures()
  failures$horizon <- numeric()
  record_failure <- function(error, stage, look, n, horizon) {
    if (!is.null(error)) {
      failures <<- rbind(
        failures,
        data.frame(
          trial = trial,
          look = as.integer(look),
          planned_N = as.integer(n),
          stage = stage,
          error_class = error$error_class,
          message = error$message,
          horizon = horizon
        )
      )
    }
  }
  generation <- args
  # Convert once from the fixed dropout reference, never from a look's endpoint.
  if (args$followup_horizon != args$dropout_horizon) {
    generation$prop_loss <- -expm1(
      log1p(-args$prop_loss) * args$followup_horizon / args$dropout_horizon
    )
  }
  generated <- capture_path_stage(
    stream,
    call_path_function(sim_comp_data, generation)
  )
  record_failure(
    generated$error,
    "generation",
    NA_integer_,
    args$N_total,
    args$followup_horizon
  )
  interims <- vector("list", length(looks))
  finals <- list()
  if (is.null(generated$value)) {
    return(list(
      trial = trial,
      interims = interims,
      finals = NULL,
      failures = failures,
      generation_warnings = generated$warnings
    ))
  }
  data <- generated$value
  stopped <- FALSE
  stop_index <- length(cohorts)
  for (i in seq_along(cohorts)) {
    n <- cohorts[i]
    if (i <= length(looks)) {
      evidence <- vector("list", length(horizons))
      needed <- if (is.null(schedule)) {
        seq_along(horizons)
      } else {
        unique(c(match(schedule$horizon[i], horizons), final_index))
      }
      for (j in seq_along(horizons)) {
        # Reserve every slot, including unused calculations in direct simulation.
        stream <- parallel::nextRNGSubStream(stream)
        if (stopped || !j %in% needed) {
          next
        }
        stage <- capture_path_stage(
          stream,
          horizon_interim_evidence(data, configs[[j]], i, n)
        )
        record_failure(stage$error, "interim", i, n, horizons[j])
        if (!is.null(stage$value)) {
          for (cohort in c("current", "maximum")) {
            record_failure(
              attr(stage$value$scores[[cohort]], "failure", exact = TRUE),
              paste0("interim_", cohort),
              i,
              n,
              horizons[j]
            )
          }
          stage$value$diagnostics$warnings <- unique(c(
            stage$value$diagnostics$warnings,
            stage$warnings
          ))
          stage$value$diagnostics$prior <- NULL
        }
        evidence[j] <- list(stage$value)
      }
      interims[[i]] <- evidence
      if (!is.null(schedule) && !stopped) {
        combined <- combine_horizon_evidence(
          evidence[[match(schedule$horizon[i], horizons)]],
          evidence[[final_index]],
          schedule$horizon[i],
          args$final_horizon
        )
        decision <- tryCatch(
          horizon_interim_decision(combined, rules, i),
          error = function(e) "failed"
        )
        if (decision != "continue") {
          stopped <- TRUE
          stop_index <- i
        }
      }
    }
    for (j in seq_along(horizons)) {
      stream <- parallel::nextRNGSubStream(stream)
      cohort <- censor_horizon_data(
        data[data$id <= n, , drop = FALSE],
        horizons[j]
      )
      cohort$subject_impute_success <- cohort$event == 0 &
        cohort$time < horizons[j]
      needed <- is.null(schedule) ||
        (i == stop_index &&
          horizons[j] ==
            if (i <= length(looks)) {
              schedule$horizon[i]
            } else {
              args$final_horizon
            }) ||
        (i == length(cohorts) && j == final_index)
      stage <- if (needed) {
        capture_path_stage(
          stream,
          call_path_function(analyse_final, configs[[j]], data_in = cohort)
        )
      } else {
        list(value = NULL, error = NULL, warnings = character())
      }
      if (
        !is.null(stage$value) &&
          (length(stage$value) != 2L ||
            !is.finite(stage$value[1]) ||
            stage$value[1] < 0 ||
            stage$value[1] > 1)
      ) {
        stage$error <- list(
          error_class = "goldilocks_nonfinite_final",
          message = "Final analysis did not return a finite probability"
        )
        stage$value <- NULL
      }
      record_failure(stage$error, "final", i, n, horizons[j])
      final <- if (is.null(stage$value)) c(NA_real_, NA_real_) else stage$value
      collection <- trial_calendar_metrics(
        data[data$id <= n, , drop = FALSE],
        args$followup_horizon
      )
      finals[[length(finals) + 1L]] <- data.frame(
        look = as.integer(i),
        planned_N = as.integer(n),
        horizon = horizons[j],
        N_treatment = sum(cohort$treatment == 1),
        N_control = sum(cohort$treatment == 0),
        post_prob_ha = final[1],
        est_final = final[2],
        trial_calendar_metrics(cohort, horizons[j]),
        collection_ready_time = collection$analysis_ready_time,
        planned_collection_time = collection$planned_completion_time,
        collection_person_time = collection$followup_person_time,
        warning_messages = paste(stage$warnings, collapse = " | ")
      )
    }
  }
  list(
    trial = trial,
    interims = interims,
    finals = dplyr::bind_rows(finals),
    failures = failures,
    generation_warnings = generated$warnings
  )
}

#' Print trial paths with multiple endpoint horizons
#' @param x A `goldilocks_horizon_paths` object.
#' @param ... Unused arguments.
#' @return The object, invisibly.
#' @export
print.goldilocks_horizon_paths <- function(x, ...) {
  args <- validate_horizon_arguments(attr(x, "arguments", exact = TRUE))
  cat("Goldilocks endpoint-horizon paths\n")
  cat(
    length(x$trials),
    "trials;",
    length(args$interim_look),
    "interim looks; horizons:",
    paste(args$endpoint_horizons, collapse = ", "),
    "\n"
  )
  cat("Maximum N:", args$N_total, "; final endpoint:", args$final_horizon, "\n")
  cat("Failed calculations:", nrow(x$failures), "\n")
  cat("Memory:", format(utils::object.size(x), units = "auto"), "\n")
  invisible(x)
}
