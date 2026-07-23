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
#'   processing. Defaults to `1L` (serial execution).
#' @param backend character. Parallel backend. "auto" (the default) uses
#'   serial execution for `ncores = 1`, the existing fork backend on Unix-like
#'   platforms, and a PSOCK cluster on Windows. "fork", "psock", and
#'   "sequential" select a backend explicitly.
#' @param seed optional integer. Seed used to generate independent per-trial
#'   `"L'Ecuyer-CMRG"` random-number streams. The default, `NULL`,
#'   does not reset the global RNG state, preserving the usual unseeded
#'   simulation behavior.
#' @param return_trace logical. Should the compact interim decision trace from
#'   every simulated trial be retained? The default, `FALSE`, preserves the
#'   compact output. When `TRUE`, the returned list also contains a `traces`
#'   data frame with a `trial` column linking each trace row to the corresponding
#'   row of `sims`.
#'
#' @details This is basically a wrapper function for
#'   [survival_adapt()], whereby we repeatedly run the function for independent
#'   trials (all with the same input design parameters and treatment effect).
#'
#'   To use multiple cores (where available), the argument `ncores`
#'   can be increased from the default of 1. The default `backend = "auto"`
#'   uses [pbmcapply::pbmclapply()] on Unix-like platforms and a PSOCK cluster
#'   on Windows, where forked processes are unavailable. Set `backend`
#'   explicitly to compare backends or to require serial execution.
#'
#'   Set `seed` to make `sim_trials()` reproducible. When a seed is
#'   supplied, `sim_trials()` first generates one independent
#'   `"L'Ecuyer-CMRG"` stream for each simulated trial, then each call to
#'   [survival_adapt()] runs with its own per-trial stream. This avoids
#'   reusing the same random-number stream across workers when
#'   `ncores > 1`, and produces identical seeded results across supported
#'   backends. A seeded call restores the caller's RNG state on exit. With
#'   `seed = NULL`, the function uses and advances R's current global RNG
#'   state.
#'
#' @return A list containing `sims`, a data frame with one row per simulated
#'   trial, and `call`. When `return_trace = TRUE`, the list also contains
#'   `traces`, a data frame with one row per completed interim look and a
#'   `trial` identifier. See [survival_adapt()] for details of the summary and
#'   trace columns.
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
#'   prior = c(0.1, 0.1),
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
  prior = c(0.1, 0.1),
  bin_prior = c(1, 1),
  bin_method = "mc",
  block = 2,
  rand_ratio = c(1, 1),
  prop_loss = 0,
  alternative = "greater",
  h0 = 0,
  Fn = 0.05,
  Sn = 0.9,
  prob_ha = 0.95,
  N_impute = 10,
  N_mcmc = 10,
  N_trials = 10,
  method = "logrank",
  imputed_final = FALSE,
  empty_interval = c("propagate", "prior", "error"),
  return_trace = FALSE,
  ncores = 1L,
  backend = c("auto", "fork", "psock", "sequential"),
  seed = NULL,
  binary_imputation = c("event-time", "bernoulli")
) {
  Call <- match.call()
  empty_interval <- match.arg(empty_interval)
  binary_imputation <- match.arg(binary_imputation)
  backend <- match.arg(backend)

  validate_positive_integer_scalar(N_trials, "N_trials")
  validate_logical_scalar(return_trace, "return_trace")

  validate_positive_integer_scalar(ncores, "ncores")
  backend <- resolve_sim_backend(backend, ncores)

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

    trial_streams <- make_rng_streams(seed, N_trials)
  } else {
    trial_streams <- NULL
  }

  survival_adapt_fn <- if (backend == "psock") {
    make_psock_callable("survival_adapt")
  } else {
    survival_adapt
  }

  survival_adapt_wrapper <- function(x) {
    if (!is.null(trial_streams)) {
      assign(".Random.seed", trial_streams[[x]], envir = .GlobalEnv)
    }
    survival_adapt_fn(
      hazard_treatment = hazard_treatment,
      hazard_control = hazard_control,
      cutpoints = cutpoints,
      N_total = N_total,
      lambda = lambda,
      lambda_time = lambda_time,
      interim_look = interim_look,
      end_of_study = end_of_study,
      prior = prior,
      bin_prior = bin_prior,
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
      method = method,
      imputed_final = imputed_final,
      empty_interval = empty_interval,
      return_trace = return_trace
    )
  }

  trial_index <- seq_len(N_trials)
  trial_results <- switch(
    backend,
    sequential = lapply(trial_index, survival_adapt_wrapper),
    fork = pbmclapply(trial_index, survival_adapt_wrapper, mc.cores = ncores),
    psock = {
      cluster <- parallel::makeCluster(ncores)
      on.exit(parallel::stopCluster(cluster), add = TRUE)
      parallel::parLapply(cluster, trial_index, survival_adapt_wrapper)
    }
  )

  if (return_trace) {
    sims <- bind_rows(lapply(trial_results, function(x) x$summary))
    traces <- bind_rows(lapply(seq_along(trial_results), function(i) {
      trace <- trial_results[[i]]$trace
      trace$trial <- rep.int(i, nrow(trace))
      trace[c("trial", setdiff(names(trace), "trial"))]
    }))
    out <- list(sims = sims, traces = traces, call = Call)
  } else {
    sims <- bind_rows(trial_results)
    out <- list(sims = sims, call = Call)
  }

  return(out)
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
#'   provide imported functions and compiled routines.
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
