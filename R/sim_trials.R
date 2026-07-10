#' @title Simulate one or more clinical trials subject to known design
#'   parameters and treatment effect
#'
#' @description Simulate multiple clinical trials with fixed input parameters,
#'   and tidily extract the relevant data to generate operating characteristics.
#'
#' @inheritParams survival_adapt
#' @inheritParams sim_comp_data
#' @param N_trials integer. Number of trials to simulate.
#' @param ncores integer. Number of cores to use for parallel processing.
#' @param seed optional integer. Seed used to generate independent per-trial
#'   `"L'Ecuyer-CMRG"` random-number streams. The default, `NULL`,
#'   does not reset the global RNG state, preserving the usual unseeded
#'   simulation behavior.
#'
#' @details This is basically a wrapper function for
#'   [survival_adapt()], whereby we repeatedly run the function for independent
#'   trials (all with the same input design parameters and treatment effect).
#'
#'   To use multiple cores (where available), the argument `ncores`
#'   can be increased from the default of 1. Note: on Windows machines, it is
#'   not possible to use [parallel::mclapply()] with `ncores` \eqn{> 1}.
#'
#'   Set `seed` to make `sim_trials()` reproducible. When a seed is
#'   supplied, `sim_trials()` first generates one independent
#'   `"L'Ecuyer-CMRG"` stream for each simulated trial, then each call to
#'   [survival_adapt()] runs with its own per-trial stream. This avoids
#'   reusing the same random-number stream across workers when
#'   `ncores > 1`. With `seed = NULL`, the function uses R's current
#'   global RNG state.
#'
#' @return Data frame with 1 row per simulated trial and columns for key summary
#'   statistics. See [survival_adapt()] for details of what is returned in each
#'   row.
#'
#' @importFrom parallel detectCores
#' @importFrom pbmcapply pbmclapply
#' @export
#'
#' @examples
#' hc <- prop_to_haz(c(0.20, 0.30), c(0, 12), 36)
#' ht <- prop_to_haz(c(0.05, 0.15), c(0, 12), 36)
#'
#' out <- sim_trials(
#'   hazard_treatment = ht,
#'   hazard_control = hc,
#'   cutpoints = c(0, 12),
#'   N_total = 600,
#'   lambda = 20,
#'   lambda_time = 0,
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
#'   seed = 123)

sim_trials <- function(
  hazard_treatment,
  hazard_control = NULL,
  cutpoints = 0,
  N_total,
  lambda = 0.3,
  lambda_time = 0,
  interim_look = NULL,
  end_of_study,
  prior = c(0.1, 0.1),
  bin_prior = c(1, 1),
  bin_method = "mc",
  bin_N = 10000,
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
  ncores = 1L,
  seed = NULL
) {
  Call <- match.call()
  empty_interval <- match.arg(empty_interval)

  validate_positive_integer_scalar(N_trials, "N_trials")

  # Check: missing 'ncores' defaults to maximum available (spare 1)
  if (missing(ncores)) {
    ncores <- default_ncores()
  } else {
    validate_positive_integer_scalar(ncores, "ncores")
  }

  # Check: cannot specify <1 core
  if (ncores < 1) {
    warning("Must use at least 1 core... setting ncores = 1")
    ncores <- 1L
  }

  # Check: if Windows and if ncores = 1
  if (.Platform$OS.type == "Windows" & ncores > 1L) {
    message("On Windows machines it is required that ncores = 1L")
    ncores <- 1
  }

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
    trial_streams <- make_rng_streams(seed, N_trials)
  } else {
    trial_streams <- NULL
  }

  survival_adapt_wrapper <- function(x) {
    if (!is.null(trial_streams)) {
      assign(".Random.seed", trial_streams[[x]], envir = .GlobalEnv)
    }
    survival_adapt(
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
      bin_N = bin_N,
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
      empty_interval = empty_interval
    )
  }

  sims <- pbmclapply(1:N_trials, survival_adapt_wrapper, mc.cores = ncores)

  sims <- bind_rows(sims)
  out <- list(sims = sims, call = Call)

  return(out)
}

#' @title Determine a default core count
#'
#' @description Reserves one detected logical core for the system and falls
#'   back to serial execution when the core count is unavailable.
#'
#' @noRd
default_ncores <- function(
  detected = parallel::detectCores(logical = TRUE)
) {
  if (
    length(detected) != 1 ||
      !is.numeric(detected) ||
      is.na(detected) ||
      !is.finite(detected) ||
      detected < 2
  ) {
    return(1L)
  }

  max(1L, as.integer(detected) - 1L)
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
