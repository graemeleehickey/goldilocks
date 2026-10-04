#' Simulate complete trial paths for repeated threshold calibration
#'
#' @description Generates one reusable simulation bank under fixed design and
#'   data-generating assumptions. Every planned interim look is evaluated,
#'   regardless of hypothetical stopping decisions. Use [apply_stopping_rules()]
#'   to calibrate `Fn`, `Sn`, `Qn`, and `prob_ha` without new simulations.
#'
#' @inheritParams sim_trials
#' @param trial_offset A non-negative integer. Trial identifiers start at
#'   `trial_offset + 1`. With an explicit seed, consecutive offsets allow a large
#'   bank to be generated, saved, and evaluated in separate batches using the
#'   same streams as a single call. Keep all other simulation settings fixed.
#'
#' @details The bank retains the success score (posterior probability or
#'   frequentist `1 - p`) of every predictive replicate, for the current and
#'   maximum cohorts, plus posterior Monte Carlo counts when applicable. It also
#'   retains a final analysis and calendar metrics for each possible stopping
#'   cohort. Final analyses use eventual observed follow-up, including the
#'   specified dropout and final-imputation policy. Interim evidence only uses
#'   information available at its calendar cut.
#'
#'   The data-generating model, priors, method, null margin, endpoint horizon,
#'   allocation, look schedule, maximum sample size, `N_impute`, and `N_mcmc`
#'   are fixed within a bank. Changing them requires a new bank. Interim looks
#'   must be strictly below `N_total`; the maximum-N analysis is after follow-up.
#'
#'   A bank costs more than one early-stopping simulation but can be reused for
#'   many candidate rules. Score storage alone requires approximately
#'   `16 * N_trials * length(interim_look) * N_impute` bytes, or 24 bytes per
#'   trial/look/imputation for Monte Carlo Bayesian methods, plus object and
#'   diagnostic overhead. Large estimates are reported before simulation. Banks
#'   are held in memory; use `trial_offset` and [saveRDS()] for separate batches.
#'
#'   Each trial and stage receives a separate L'Ecuyer-CMRG stream or substream.
#'   An explicit seed preserves the caller's RNG state and gives identical paths
#'   across supported backends. With `seed = NULL`, exactly one integer seed is
#'   drawn from the caller's RNG; all backends use that seed. Seeded bank results
#'   need not equal [sim_trials()] results, whose draw scheduling is unchanged.
#'
#'   Failures are retained at the calculation level. A failed predictive
#'   replicate invalidates its entire current- or maximum-cohort probability;
#'   replicates are never silently dropped. Later looks and cohort final
#'   analyses are still attempted. Rule application determines which failures
#'   prevent a particular candidate from being evaluated.
#'
#' @return A `goldilocks_paths` object containing `trials` (one record per
#'   requested trial), `failures`, `schema_version`, and `call`. Each trial has
#'   `interims` with unthresholded `scores`, observed `context`, and diagnostics,
#'   and `finals` with one row per possible stopping cohort. Attributes retain
#'   evaluated `arguments`, `prior_design`, `enrollment_design`, `rng_metadata`,
#'   `parallel_metadata`, `runtime_metadata`, and `score_storage_bytes`. Save the whole object with
#'   [saveRDS()] to retain its metadata and posterior counts.
#'
#' @seealso [apply_stopping_rules()], [summarise_sims()]
#' @export
#' @examples
#' paths <- sim_trial_paths(
#'   hazard_treatment = 0.05, hazard_control = 0.08,
#'   N_total = 60, interim_look = c(20, 40), lambda = 5,
#'   end_of_study = 12, alternative = "less", method = "bayes-bin",
#'   bin_method = "quadrature", N_impute = 5, N_trials = 2, seed = 123
#' )
#' result <- apply_stopping_rules(paths, Fn = 0.05, Sn = 0.9, prob_ha = 0.975)
#' summarise_sims(result)
sim_trial_paths <- function(
  hazard_treatment,
  hazard_control = NULL,
  cutpoints = NULL,
  N_total,
  lambda = 0.3,
  lambda_time = NULL,
  interim_look = NULL,
  end_of_study,
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
  binary_imputation = c("event-time", "bernoulli"),
  prior_surv_final = prior_surv,
  generation_cutpoints = cutpoints,
  rmst_tau = end_of_study,
  trial_offset = 0L
) {
  Call <- match.call()
  args <- capture_arguments(sim_trial_paths, environment())
  args$backend <- match.arg(backend)
  args$empty_interval <- match.arg(empty_interval)
  args$binary_imputation <- match.arg(binary_imputation)
  args <- validate_path_arguments(args)
  config <- args
  config$single_arm <- is.null(args$hazard_control)
  if (!args$method %in% c("bayes-surv", "bayes-bin")) {
    config$N_mcmc <- 1L
  }
  execution <- resolve_sim_execution(args$backend, args$ncores, args$N_trials)
  inner_mc <- args$method == "bayes-surv" ||
    (args$method == "bayes-bin" && args$bin_method == "mc")
  score_bytes <- as.double(N_trials) *
    length(interim_look) *
    N_impute *
    if (inner_mc) 24 else 16
  if (score_bytes >= 512 * 1024^2) {
    message(sprintf(
      "Predictive scores require approximately %.1f MiB before diagnostics and object overhead; consider separate batches with 'trial_offset'.",
      score_bytes / 1024^2
    ))
  }

  caller_kind <- RNGkind()
  stream_seed <- if (is.null(seed)) {
    sample.int(.Machine$integer.max - 1L, 1L)
  } else {
    seed
  }
  restore_rng <- preserve_path_rng()
  on.exit(restore_rng(), add = TRUE)
  trial_ids <- as.integer(trial_offset + seq_len(N_trials))
  streams <- make_rng_streams(stream_seed, max(trial_ids))[trial_ids]
  simulate <- if (execution$backend == "psock") {
    make_psock_callable("simulate_trial_path")
  } else {
    simulate_trial_path
  }
  worker <- function(i) simulate(config, streams[[i]], trial_ids[i])
  indices <- seq_len(N_trials)
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
  # The bank holds the common failure table once, rather than per trial as well.
  trials <- lapply(trials, function(x) {
    x$failures <- NULL
    x
  })
  out <- structure(
    list(
      trials = trials,
      failures = failures,
      schema_version = 1L,
      call = Call
    ),
    class = "goldilocks_paths"
  )
  attr(out, "arguments") <- args
  attr(out, "enrollment_design") <- call_path_function(
    new_enrollment_design,
    config
  )
  attr(out, "prior_design") <- dplyr::bind_rows(
    call_path_function(gamma_prior_diagnostics, config, stage = "interim"),
    call_path_function(
      gamma_prior_diagnostics,
      config,
      prior_surv = config$prior_surv_final,
      stage = "final"
    )
  )
  attr(out, "score_storage_bytes") <- score_bytes
  attr(out, "runtime_metadata") <- list(
    package_version = as.character(utils::packageVersion("goldilocks")),
    R_version = as.character(getRversion())
  )
  attr(out, "rng_metadata") <- list(
    caller_kind = caller_kind,
    stream_kind = "L'Ecuyer-CMRG",
    seed_policy = if (is.null(seed)) {
      "caller_derived_bank"
    } else {
      "explicit_preserve_caller"
    },
    stream_seed = stream_seed,
    stream_layout = "trial-stage-v1",
    backend = execution$backend,
    ncores = as.integer(execution$workers)
  )
  attr(out, "parallel_metadata") <- list(
    requested_backend = args$backend,
    backend = execution$backend,
    selection_reason = execution$reason,
    requested_ncores = as.integer(ncores),
    workers = as.integer(execution$workers),
    tasks = as.integer(N_trials)
  )
  if (nrow(failures)) {
    warning(
      nrow(failures),
      " path calculations failed; retained in 'paths$failures'. ",
      "Candidate rules determine which trials are evaluable.",
      call. = FALSE
    )
  }
  out
}

# Validate fixed settings before workers or random-number generation begin.
validate_path_arguments <- function(args) {
  args$method <- normalize_analysis_method(args$method)
  for (name in c("N_total", "N_trials", "N_impute", "N_mcmc", "ncores")) {
    validate_positive_integer_scalar(args[[name]], name)
  }
  validate_nonnegative_integer_scalar(args$trial_offset, "trial_offset")
  if (as.double(args$trial_offset) + args$N_trials > .Machine$integer.max) {
    stop(
      "'trial_offset + N_trials' exceeds the supported trial identifier range"
    )
  }
  if (!is.null(args$seed)) {
    validate_nonnegative_integer_scalar(args$seed, "seed")
    if (args$seed > .Machine$integer.max) {
      stop("'seed' exceeds .Machine$integer.max")
    }
  }
  single_arm <- is.null(args$hazard_control)
  if (!single_arm) {
    args$rand_ratio <- validate_randomization_args(
      args$N_total,
      args$block,
      args$rand_ratio,
      allocation_name = "rand_ratio"
    )
  }
  args$prop_loss <- normalize_prop_loss(args$prop_loss, single_arm)
  validate_single_probability(
    args$mc_conf_level,
    "mc_conf_level",
    upper_open = TRUE
  )
  if (args$mc_conf_level <= 0.5) {
    stop("'mc_conf_level' must be greater than 0.5 and less than 1")
  }
  validate_cutpoints(args$cutpoints)
  validate_endpoint_time(args$end_of_study, args$cutpoints, "end_of_study")
  validate_cutpoints(args$generation_cutpoints, "generation_cutpoints")
  validate_endpoint_time(
    args$end_of_study,
    args$generation_cutpoints,
    "end_of_study",
    "generation_cutpoints"
  )
  validate_piecewise_hazard(
    args$hazard_treatment,
    args$generation_cutpoints,
    "hazard_treatment",
    "generation_cutpoints"
  )
  if (!single_arm) {
    validate_piecewise_hazard(
      args$hazard_control,
      args$generation_cutpoints,
      "hazard_control",
      "generation_cutpoints"
    )
  }
  validate_enrollment_schedule(args$lambda, args$lambda_time, args$N_total)
  for (name in c("prior_surv", "prior_surv_final")) {
    args[[name]] <- normalize_gamma_prior(
      args[[name]],
      length(args$cutpoints) + 1L,
      single_arm,
      name
    )
  }
  validate_analysis_configuration(
    args$method,
    args$alternative,
    single_arm,
    args$imputed_final
  )
  validate_final_imputation(
    args$method,
    args$imputed_final,
    has_missing_outcomes = any(args$prop_loss > 0),
    N_impute = args$N_impute
  )
  if (!is.null(args$interim_look)) {
    validate_interim_looks(
      args$interim_look,
      args$N_total,
      if (single_arm) NULL else max(args$block)
    )
  }
  validate_h0(args$h0, args$method, single_arm)
  if (args$method == "rmst") {
    validate_rmst_args(args$rmst_tau, args$end_of_study, args$h0)
  }
  if (args$method == "bayes-bin") {
    validate_bayes_binomial_args(args$prior_bin, args$bin_method, args$N_mcmc)
  }
  args
}

# Forward only arguments declared by the shared calculation, with explicit overrides.
call_path_function <- function(fun, config, ...) {
  args <- config[intersect(names(formals(fun)), names(config))]
  extra <- list(...)
  args[names(extra)] <- extra
  do.call(fun, args)
}

preserve_path_rng <- function() {
  kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  saved <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  function() {
    do.call(RNGkind, as.list(kind))
    if (had_seed) {
      assign(".Random.seed", saved, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }
}

capture_path_stage <- function(stream, expr) {
  assign(".Random.seed", stream, envir = .GlobalEnv)
  warnings <- character()
  value <- tryCatch(
    withCallingHandlers(expr, warning = function(w) {
      warnings <<- unique(c(warnings, conditionMessage(w)))
      invokeRestart("muffleWarning")
    }),
    error = identity
  )
  list(
    value = if (inherits(value, "error")) NULL else value,
    error = if (inherits(value, "error")) {
      list(error_class = class(value)[1L], message = conditionMessage(value))
    } else {
      NULL
    },
    warnings = warnings
  )
}

empty_path_failures <- function() {
  data.frame(
    trial = integer(),
    look = integer(),
    planned_N = integer(),
    stage = character(),
    error_class = character(),
    message = character()
  )
}

simulate_trial_path <- function(config, stream, trial) {
  looks <- config$interim_look
  cohorts <- c(looks, config$N_total)
  failures <- empty_path_failures()
  record_failure <- function(error, stage, look, n) {
    if (!is.null(error)) {
      failures <<- rbind(
        failures,
        data.frame(
          trial = trial,
          look = as.integer(look),
          planned_N = as.integer(n),
          stage = stage,
          error_class = error$error_class,
          message = error$message
        )
      )
    }
  }
  generated <- capture_path_stage(
    stream,
    call_path_function(sim_comp_data, config)
  )
  record_failure(generated$error, "generation", NA_integer_, config$N_total)
  if (is.null(generated$value)) {
    return(list(
      trial = trial,
      interims = vector("list", length(looks)),
      finals = NULL,
      failures = failures,
      generation_warnings = generated$warnings
    ))
  }
  data_total <- generated$value
  interims <- vector("list", length(looks))
  finals <- vector("list", length(cohorts))
  for (i in seq_along(cohorts)) {
    n <- cohorts[i]
    if (i <= length(looks)) {
      stream <- parallel::nextRNGSubStream(stream)
      stage <- capture_path_stage(
        stream,
        call_path_function(
          evaluate_interim_evidence,
          config,
          data_interim = prepare_simulated_interim(
            data_total,
            n,
            config$end_of_study
          ),
          look = i,
          planned_N = n,
          calendar_time = data_total$enrollment[n],
          active_followup = active_followup_at(
            data_total,
            data_total$enrollment[n]
          ),
          check_futility = TRUE,
          isolate_errors = TRUE
        )
      )
      record_failure(stage$error, "interim", i, n)
      if (!is.null(stage$value)) {
        for (cohort in c("current", "maximum")) {
          record_failure(
            attr(stage$value$scores[[cohort]], "failure", exact = TRUE),
            paste0("interim_", cohort),
            i,
            n
          )
        }
        stage$value$diagnostics$warnings <- unique(c(
          stage$value$diagnostics$warnings,
          stage$warnings
        ))
        # Resolved priors are stored once in the bank metadata.
        stage$value$diagnostics$prior <- NULL
      }
      interims[i] <- list(stage$value)
    }
    data_final <- data_total[data_total$id <= n, , drop = FALSE]
    data_final$subject_impute_success <- data_final$event == 0 &
      data_final$time < config$end_of_study
    stream <- parallel::nextRNGSubStream(stream)
    stage <- capture_path_stage(
      stream,
      call_path_function(analyse_final, config, data_in = data_final)
    )
    if (
      !is.null(stage$value) &&
        (length(stage$value) != 2L ||
          !is.finite(stage$value[1]) ||
          stage$value[1] < 0 ||
          stage$value[1] > 1)
    ) {
      stage$error <- list(
        error_class = "goldilocks_nonfinite_final",
        message = "Final analysis did not return finite results"
      )
      stage$value <- NULL
    }
    record_failure(stage$error, "final", i, n)
    final <- if (is.null(stage$value)) c(NA_real_, NA_real_) else stage$value
    finals[[i]] <- data.frame(
      look = as.integer(i),
      planned_N = as.integer(n),
      N_treatment = sum(data_final$treatment == 1),
      N_control = sum(data_final$treatment == 0),
      post_prob_ha = final[1],
      est_final = final[2],
      trial_calendar_metrics(data_final, config$end_of_study),
      warning_messages = paste(stage$warnings, collapse = " | ")
    )
  }
  list(
    trial = trial,
    interims = interims,
    finals = dplyr::bind_rows(finals),
    failures = failures,
    generation_warnings = generated$warnings
  )
}

#' Print a complete trial-path simulation bank
#' @param x A `goldilocks_paths` object.
#' @param ... Unused arguments.
#' @return The object, invisibly.
#' @export
print.goldilocks_paths <- function(x, ...) {
  args <- attr(x, "arguments", exact = TRUE)
  cat("Goldilocks trial paths\n")
  cat(
    length(x$trials),
    "trials;",
    length(args$interim_look),
    "interim looks; maximum N =",
    args$N_total,
    "\n"
  )
  cat("Failed calculations:", nrow(x$failures), "\n")
  cat("Memory:", format(utils::object.size(x), units = "auto"), "\n")
  invisible(x)
}
