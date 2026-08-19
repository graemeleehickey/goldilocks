#' @title Simulate one or more clinical trials subject to known design
#'   parameters and treatment effect
#'
#' @description Simulate multiple clinical trials with fixed input parameters,
#'   and tidily extract the relevant data to generate operating characteristics.
#'
#' @inheritParams survival_adapt
#' @inheritParams sim_comp_data
#' @param N_trials integer. Number of trials to simulate.
#' @param ncores positive integer. Number of cores to use for parallel
#'   processing. Defaults to `1L` (serial execution). Package-owned worker
#'   pools are capped at the number of trials. With `backend = "auto"`, at
#'   least two trials per worker are required so that small workloads avoid
#'   parallel startup overhead.
#' @param backend character. Parallel backend. "auto" (the default) uses
#'   serial execution for `ncores = 1` or fewer than four trials. Otherwise it
#'   uses the existing fork backend on Unix-like platforms and a PSOCK cluster
#'   on Windows, with at least two trials assigned per worker. "fork",
#'   "psock", and "sequential" select a backend explicitly and bypass this
#'   workload crossover.
#' @param seed optional integer. Seed used to generate independent per-trial
#'   `"L'Ecuyer-CMRG"` random-number streams. The default, `NULL`,
#'   does not reset the global RNG state, preserving the usual unseeded
#'   simulation behavior.
#' @param return_trace logical. Should the compact interim decision trace from
#'   every simulated trial be retained? The default, `FALSE`, preserves the
#'   compact output. When `TRUE`, the returned list also contains a `traces`
#'   data frame with a `trial` column linking each trace row to the corresponding
#'   original simulated trial.
#'
#' @details This is basically a wrapper function for
#'   [survival_adapt()], whereby we repeatedly run the function for independent
#'   trials (all with the same input design parameters and treatment effect).
#'
#'   To use multiple cores (where available), the argument `ncores`
#'   can be increased from the default of 1. The default `backend = "auto"`
#'   stays sequential for fewer than four trials and otherwise uses at most one
#'   worker per two trials. This simple workload heuristic amortizes process
#'   startup without trying to predict the cost of an individual trial. When
#'   it selects parallel execution, it uses [pbmcapply::pbmclapply()] on
#'   Unix-like platforms and a PSOCK cluster on Windows, where forked processes
#'   are unavailable. Set `backend` explicitly to bypass the crossover.
#'
#'   A PSOCK cluster created by `sim_trials()` is always stopped before the
#'   function returns, including when worker evaluation fails. PSOCK tasks use
#'   static worker chunks, so the invariant simulation design is serialized
#'   once per active worker rather than once per trial.
#'
#'   Errors raised by an individual [survival_adapt()] call are isolated so
#'   other trials can finish. Failed trials are excluded from `sims`, recorded
#'   in `failures` with their trial number, error class, and message, and
#'   reported together in one warning. If every requested trial fails,
#'   `sim_trials()` stops and attaches the same failure table to the error as
#'   `failures`. With a supplied `seed`, the original call and failed trial
#'   number reproduce the same per-trial random-number stream.
#'
#'   Set `seed` to make `sim_trials()` reproducible. When a seed is
#'   supplied, `sim_trials()` first generates one independent
#'   `"L'Ecuyer-CMRG"` stream for each simulated trial, then each call to
#'   [survival_adapt()] runs with its own per-trial stream. This avoids
#'   reusing the same random-number stream across workers when
#'   `ncores > 1`, and produces identical seeded results across supported
#'   backends. A seeded call restores the caller's RNG state on exit. With
#'   `seed = NULL`, sequential execution uses and advances R's current global
#'   RNG state directly. Fork execution delegates stream setup to
#'   [pbmcapply::pbmclapply()]. PSOCK execution draws one integer seed from the
#'   caller's current RNG state, thereby advancing that state predictably, and
#'   expands it into independent per-trial `"L'Ecuyer-CMRG"` streams. Resetting
#'   the caller to the same state therefore reproduces an unseeded PSOCK call,
#'   while consecutive calls do not reuse streams.
#'
#' @return A list containing `sims`, a data frame with one row per successfully
#'   simulated trial; `failures`, a data frame with columns `trial`,
#'   `error_class`, and `message`; and `call`. When `return_trace = TRUE`, the
#'   list also contains `traces`, a data frame with one row per completed
#'   interim look and a `trial` identifier. Per-trial calendar-time metrics are
#'   always retained in `sims`; traces additionally retain calendar time and
#'   active follow-up at each look. See [survival_adapt()] for details of the
#'   summary and trace columns, and [summarise_calendar_time()] for wide
#'   operating-characteristic tables. The returned object also retains the
#'   normalized `decision_design` attribute from [survival_adapt()]. An
#'   `rng_metadata` attribute records the caller RNG kind, effective backend,
#'   seed policy, and stream seed when applicable. A deterministic
#'   `parallel_metadata` attribute records the requested and effective backend
#'   and worker counts, task count, and backend-selection reason. An `arguments`
#'   attribute contains a named list of all evaluated argument values, including
#'   defaults. It can be saved with [saveRDS()] and supplied to a later call
#'   with `do.call(sim_trials, attr(result, "arguments"))`.
#'
#' @importFrom pbmcapply pbmclapply
#' @export
#'
#' @examples
#' hc <- prop_to_haz(c(0.20, 0.30), 12, 36)
#' ht <- prop_to_haz(c(0.05, 0.15), 12, 36)
#'
#' out <- sim_trials(
#'   hazard_treatment = ht,
#'   hazard_control = hc,
#'   cutpoints = 12,
#'   N_total = 600,
#'   lambda = 20,
#'   lambda_time = NULL,
#'   interim_look = c(400, 500),
#'   end_of_study = 36,
#'   prior_surv = c(0.1, 0.1),
#'   block = 2,
#'   rand_ratio = c(1, 1),
#'   prop_loss = 0.30,
#'   alternative = "two.sided",
#'   h0 = 0,
#'   Fn = 0.05,
#'   Sn = 0.9,
#'   prob_ha = 0.975,
#'   N_impute = 5,
#'   N_mcmc = 5,
#'   method = "logrank",
#'   N_trials = 2,
#'   ncores = 1,
#'   backend = "auto",
#'   seed = 123)

sim_trials <- function(
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
  rand_ratio = c(1, 1),
  prop_loss = 0,
  alternative = "greater",
  h0 = 0,
  Fn = 0.05,
  Sn = 0.9,
  prob_ha = 0.95,
  N_impute = 500,
  N_mcmc = 1000,
  mc_conf_level = 0.95,
  N_trials = 10,
  method = "logrank",
  imputed_final = FALSE,
  empty_interval = c("prior", "propagate", "error"),
  return_trace = FALSE,
  ncores = 1L,
  backend = c("auto", "fork", "psock", "sequential"),
  seed = NULL,
  binary_imputation = c("event-time", "bernoulli"),
  prior_surv_final = prior_surv
) {
  Call <- match.call()
  Arguments <- capture_arguments(sim_trials, environment())
  empty_interval <- match.arg(empty_interval)
  binary_imputation <- match.arg(binary_imputation)
  backend <- match.arg(backend)
  Arguments$empty_interval <- empty_interval
  Arguments$binary_imputation <- binary_imputation
  Arguments$backend <- backend
  requested_backend <- backend
  caller_rng_kind <- RNGkind()

  validate_positive_integer_scalar(N_trials, "N_trials")
  validate_logical_scalar(return_trace, "return_trace")
  validate_positive_integer_scalar(N_impute, "N_impute")
  validate_positive_integer_scalar(N_mcmc, "N_mcmc")
  validate_single_probability(
    mc_conf_level,
    "mc_conf_level",
    upper_open = TRUE
  )
  if (mc_conf_level <= 0.5) {
    stop("'mc_conf_level' must be greater than 0.5 and less than 1")
  }

  validate_positive_integer_scalar(ncores, "ncores")
  execution <- resolve_sim_execution(
    backend = backend,
    ncores = ncores,
    N_trials = N_trials
  )
  backend <- execution$backend
  workers <- execution$workers

  stream_seed <- NULL
  if (!is.null(seed)) {
    if (
      length(seed) != 1 ||
        !is.numeric(seed) ||
        is.na(seed) ||
        !is.finite(seed) ||
        seed != floor(seed)
    ) {
      stop("'seed' must be NULL or a single integer value")
    }
    old_kind <- RNGkind()
    old_seed_exists <- exists(
      ".Random.seed",
      envir = .GlobalEnv,
      inherits = FALSE
    )
    if (old_seed_exists) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
    on.exit(
      {
        do.call(RNGkind, as.list(old_kind))
        if (old_seed_exists) {
          assign(".Random.seed", old_seed, envir = .GlobalEnv)
        } else if (
          exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        ) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      },
      add = TRUE
    )

    stream_seed <- seed
    trial_streams <- make_rng_streams(stream_seed, N_trials)
  } else if (backend == "psock") {
    # PSOCK workers otherwise initialize independently of the caller. Consume
    # exactly one draw from the caller's RNG, then deterministically expand it
    # into one independent stream per simulated trial.
    stream_seed <- sample.int(.Machine$integer.max - 1L, size = 1L)
    trial_streams <- make_rng_streams(stream_seed, N_trials)
  } else {
    trial_streams <- NULL
  }

  rng_metadata <- list(
    caller_kind = caller_rng_kind,
    stream_kind = if (is.null(trial_streams)) {
      caller_rng_kind[1]
    } else {
      "L'Ecuyer-CMRG"
    },
    seed_policy = if (!is.null(seed)) {
      "explicit_preserve_caller"
    } else if (backend == "psock") {
      "caller_derived_psock"
    } else {
      "caller_state"
    },
    backend = backend,
    ncores = as.integer(workers),
    requested_ncores = as.integer(ncores),
    stream_seed = stream_seed
  )

  survival_adapt_fn <- if (backend == "psock") {
    make_psock_callable("survival_adapt")
  } else {
    survival_adapt
  }

  survival_adapt_wrapper <- function(x) {
    tryCatch(
      {
        if (!is.null(trial_streams)) {
          assign(".Random.seed", trial_streams[[x]], envir = .GlobalEnv)
        }
        result <- survival_adapt_fn(
          hazard_treatment = hazard_treatment,
          hazard_control = hazard_control,
          cutpoints = cutpoints,
          N_total = N_total,
          lambda = lambda,
          lambda_time = lambda_time,
          interim_look = interim_look,
          end_of_study = end_of_study,
          prior_surv = prior_surv,
          prior_bin = prior_bin,
          bin_method = bin_method,
          binary_imputation = binary_imputation,
          block = block,
          rand_ratio = rand_ratio,
          prop_loss = prop_loss,
          alternative = alternative,
          h0 = h0,
          Fn = Fn,
          Sn = Sn,
          prob_ha = prob_ha,
          N_impute = N_impute,
          N_mcmc = N_mcmc,
          mc_conf_level = mc_conf_level,
          method = method,
          imputed_final = imputed_final,
          empty_interval = empty_interval,
          return_trace = return_trace,
          prior_surv_final = prior_surv_final
        )
        attr(result, "arguments") <- NULL
        if (inherits(result, "goldilocks_trial")) {
          attr(result$summary, "arguments") <- NULL
        }
        list(
          trial = as.integer(x),
          result = result
        )
      },
      error = function(error) {
        list(
          trial = as.integer(x),
          result = NULL,
          error_class = class(error)[1L],
          message = conditionMessage(error)
        )
      }
    )
  }

  trial_index <- seq_len(N_trials)
  trial_results <- switch(
    backend,
    sequential = lapply(trial_index, survival_adapt_wrapper),
    fork = pbmclapply(
      trial_index,
      survival_adapt_wrapper,
      mc.cores = workers
    ),
    psock = {
      active_cluster <- make_sim_cluster(workers)
      on.exit(stop_sim_cluster(active_cluster), add = TRUE)
      initialize_sim_cluster(active_cluster)
      run_sim_cluster(
        active_cluster,
        trial_index,
        survival_adapt_wrapper
      )
    }
  )

  failed <- vapply(
    trial_results,
    function(x) is.null(x$result),
    logical(1)
  )
  failed_results <- trial_results[failed]
  failures <- data.frame(
    trial = vapply(failed_results, `[[`, integer(1), "trial"),
    error_class = vapply(
      failed_results,
      `[[`,
      character(1),
      "error_class"
    ),
    message = vapply(failed_results, `[[`, character(1), "message"),
    stringsAsFactors = FALSE
  )
  if (all(failed)) {
    error <- simpleError(paste0(
      "All ",
      N_trials,
      " simulated trials failed. First error: ",
      failures$message[1L]
    ))
    error$failures <- failures
    class(error) <- c("goldilocks_all_trials_failed", class(error))
    stop(error)
  }

  successful_trials <- trial_results[!failed]
  successful_results <- lapply(successful_trials, `[[`, "result")
  if (return_trace) {
    sims <- bind_rows(lapply(successful_results, function(x) x$summary))
    traces <- bind_rows(lapply(seq_along(successful_results), function(i) {
      trace <- successful_results[[i]]$trace
      trial <- successful_trials[[i]]$trial
      trace$trial <- rep.int(trial, nrow(trace))
      trace[c("trial", setdiff(names(trace), "trial"))]
    }))
    out <- list(
      sims = sims,
      traces = traces,
      failures = failures,
      call = Call
    )
  } else {
    sims <- bind_rows(successful_results)
    out <- list(sims = sims, failures = failures, call = Call)
  }
  attr(out$sims, "arguments") <- NULL
  attr(out, "enrollment_design") <- new_enrollment_design(
    lambda = lambda,
    N_total = N_total,
    lambda_time = lambda_time,
    interim_look = interim_look,
    end_of_study = end_of_study
  )
  attr(out, "decision_design") <- attr(
    successful_results[[1]],
    "decision_design",
    exact = TRUE
  )
  attr(out, "rng_metadata") <- rng_metadata
  attr(out, "arguments") <- Arguments
  attr(out, "parallel_metadata") <- list(
    requested_backend = requested_backend,
    backend = backend,
    selection_reason = execution$reason,
    requested_ncores = as.integer(ncores),
    workers = as.integer(workers),
    tasks = as.integer(N_trials)
  )

  if (nrow(failures) > 0L) {
    warning(
      nrow(failures),
      " of ",
      N_trials,
      " simulated trials failed and were excluded. See `result$failures` ",
      "for details.",
      call. = FALSE
    )
  }

  return(out)
}

#' Resolve a trial-simulation execution plan
#'
#' @title Resolve trial-simulation execution details
#'
#' @description Applies the automatic workload crossover, caps package-owned
#'   workers at useful independent work.
#'
#' @param backend Requested backend name.
#' @param ncores Requested number of workers.
#' @param N_trials Number of independent trial tasks.
#'
#' @return A list describing the effective backend and worker pool.
#'
#' @keywords internal
#' @noRd
resolve_sim_execution <- function(backend, ncores, N_trials) {
  if (backend != "auto") {
    backend <- resolve_sim_backend(backend, ncores)
  }
  if (backend == "sequential") {
    return(list(
      backend = "sequential",
      workers = 1L,
      reason = "explicit_sequential"
    ))
  }

  if (backend == "auto") {
    auto_workers <- min(ncores, max(1L, N_trials %/% 2L))
    if (auto_workers < 2L) {
      reason <- if (ncores == 1L) {
        "single_worker"
      } else {
        "auto_small_workload"
      }
      return(list(
        backend = "sequential",
        workers = 1L,
        reason = reason
      ))
    }
    effective_backend <- resolve_sim_backend("auto", auto_workers)
    return(list(
      backend = effective_backend,
      workers = as.integer(auto_workers),
      reason = "auto_parallel"
    ))
  }

  workers <- as.integer(min(ncores, N_trials))
  list(
    backend = backend,
    workers = workers,
    reason = "explicit_parallel"
  )
}

# Small wrappers keep package-owned cluster lifecycle behavior independently
# testable without opening worker processes.
make_sim_cluster <- function(workers) {
  parallel::makeCluster(workers)
}

initialize_sim_cluster <- function(cluster) {
  source_namespace <- environment(sim_trials)
  package_name <- getNamespaceName(source_namespace)
  package_path <- getNamespaceInfo(source_namespace, "path")
  package_metadata <- file.path(package_path, "Meta", "package.rds")

  # A namespace created by pkgload::load_all() does not have an installed
  # package library to propagate. Its serialized namespace is sufficient for
  # local Unix-like PSOCK tests; installed-package execution takes this path.
  if (length(package_path) != 1L || !file.exists(package_metadata)) {
    return(invisible(cluster))
  }

  package_library <- dirname(package_path)
  worker_library_paths <- unique(c(package_library, .libPaths()))
  worker_initializer <- function(
    package_name,
    package_library,
    library_paths
  ) {
    .libPaths(library_paths)
    loadNamespace(package_name, lib.loc = package_library)
    invisible(NULL)
  }
  environment(worker_initializer) <- baseenv()

  parallel::clusterCall(
    cluster,
    worker_initializer,
    package_name,
    package_library,
    worker_library_paths
  )
  invisible(cluster)
}

run_sim_cluster <- function(cluster, trial_index, fun) {
  parallel::parLapply(cluster, trial_index, fun)
}

stop_sim_cluster <- function(cluster) {
  parallel::stopCluster(cluster)
}

#' Resolve the execution backend for trial simulation
#'
#' @title Resolve a trial-simulation backend
#'
#' @description Maps the platform-independent "auto" choice to the serial,
#'   fork, or PSOCK implementation. Forking is rejected on Windows because R
#'   does not support it there.
#'
#' @param backend Requested backend name.
#' @param ncores Number of requested workers.
#'
#' @return A single backend name.
#'
#' @keywords internal
#' @noRd
resolve_sim_backend <- function(backend, ncores) {
  if (backend == "auto") {
    if (ncores == 1L) {
      return("sequential")
    }
    if (.Platform$OS.type == "windows") {
      return("psock")
    }
    return("fork")
  }

  if (backend == "fork" && .Platform$OS.type == "windows") {
    stop("'backend = \"fork\"' is not supported on Windows")
  }

  backend
}

#' Prepare a package function for a PSOCK worker
#'
#' @title Prepare a package function for PSOCK execution
#'
#' @description Re-homes the package's R functions in a serializable
#'   environment so PSOCK workers use the same source implementation as the
#'   calling session. The package namespace remains the parent environment to
#'   provide imported functions and compiled routines. This helper is used by
#'   multi-core Windows execution, where `sim_trials()` selects PSOCK workers
#'   because forked processes are unavailable. Before the callable is sent,
#'   `initialize_sim_cluster()` gives each worker the parent package library
#'   and loads the namespace so registered compiled routines are available.
#'
#' @param name Name of the package function to prepare.
#'
#' @return A function with all package R dependencies available in its
#'   enclosing environment.
#'
#' @keywords internal
#' @noRd
make_psock_callable <- function(name) {
  source_namespace <- environment(sim_trials)
  worker_environment <- new.env(parent = source_namespace)
  function_names <- ls(source_namespace, all.names = TRUE)
  function_names <- function_names[
    vapply(
      function_names,
      function(x) is.function(get(x, envir = source_namespace)),
      logical(1)
    )
  ]

  for (function_name in function_names) {
    worker_function <- get(function_name, envir = source_namespace)
    environment(worker_function) <- worker_environment
    assign(function_name, worker_function, envir = worker_environment)
  }

  get(name, envir = worker_environment)
}

#' Generate independent random-number streams for trial simulations
#'
#' @title Create per-trial random-number streams
#'
#' @description Creates one `"L'Ecuyer-CMRG"` random-number stream per
#'   simulated trial, while preserving the caller's existing RNG kind and
#'   global `.Random.seed`. These streams are assigned inside each
#'   `survival_adapt()` call so seeded simulations are reproducible across
#'   serial and parallel execution.
#'
#' @param seed Integer seed used to initialize the stream sequence.
#' @param n Integer number of streams to generate.
#'
#' @return A list of length `n`; each element is an integer vector that can be
#'   assigned to `.Random.seed`.
#'
#' @keywords internal
#' @noRd
make_rng_streams <- function(seed, n) {
  old_kind <- RNGkind()
  old_seed_exists <- exists(
    ".Random.seed",
    envir = .GlobalEnv,
    inherits = FALSE
  )
  if (old_seed_exists) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }

  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (old_seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  })

  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)

  streams <- vector("list", n)
  streams[[1]] <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  for (i in seq_len(n - 1) + 1L) {
    streams[[i]] <- parallel::nextRNGStream(streams[[i - 1L]])
  }

  streams
}
