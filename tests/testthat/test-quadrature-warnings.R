collect_quadrature_warnings <- function(expr) {
  warnings <- list()
  value <- withCallingHandlers(
    expr,
    goldilocks_quadrature_warning = function(w) {
      warnings[[length(warnings) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = warnings)
}

quadrature_beta <- function(alpha, beta) {
  list(
    alpha = alpha,
    beta = beta,
    mean = alpha / (alpha + beta),
    variance = alpha * beta / ((alpha + beta)^2 * (alpha + beta + 1))
  )
}

quadrature_warning_args <- function() {
  # Concentrated priors reproduce the narrow-posterior integration problem
  # in small trials. Negligible hazards keep the event counts at zero.
  list(
    hazard_treatment = 1e-12,
    hazard_control = 1e-12,
    N_total = 40,
    lambda = 5,
    interim_look = c(20, 30),
    end_of_study = 12,
    prior_bin = c(1, 10001),
    prior_surv = c(1, 1e10),
    method = "bayes-bin",
    bin_method = "quadrature",
    alternative = "less",
    block = 2,
    N_impute = 3,
    N_mcmc = NULL,
    N_trials = 3,
    seed = 441,
    backend = "sequential"
  )
}

test_that("symmetry detects a known quadrature failure without replacing it", {
  set.seed(192)
  rng <- .Random.seed
  actual <- collect_quadrature_warnings(bayes_binomial_from_counts(
    events_control = 0,
    n_control = 10000,
    events_treatment = 0,
    n_treatment = 10000,
    single_arm = FALSE,
    alternative = "less",
    h0 = 0,
    prior_bin = c(1, 1),
    bin_method = "quadrature",
    N_mcmc = NULL
  ))
  # Identical Beta(1, 10001) distributions imply 0.5, but the original
  # integration misses the narrow mass and returns approximately 1.07e-8.
  expect_lt(actual$value$success, 1e-6)
  expect_equal(actual$value$effect, 0)
  expect_identical(.Random.seed, rng)
  expect_length(actual$warnings, 1)
  expect_equal(
    actual$warnings[[1]]$counts,
    c(
      flagged_calculations = 1,
      symmetry = 1,
      moments = 0
    )
  )
  expect_match(conditionMessage(actual$warnings[[1]]), "have not been changed")

  posterior <- quadrature_beta(2, 8)
  for (direction in c("less", "greater")) {
    expect_no_warning(check_beta_quadrature(
      0.5 + 5e-7,
      posterior,
      posterior,
      0,
      direction
    ))
    expect_warning(
      check_beta_quadrature(
        0.5 + 2e-6,
        posterior,
        posterior,
        0,
        direction
      ),
      class = "goldilocks_quadrature_warning"
    )
  }
})

test_that("posterior moment bounds detect errors with unequal posteriors", {
  actual <- collect_quadrature_warnings(bayes_binomial_from_counts(
    events_control = 1,
    n_control = 10000,
    events_treatment = 0,
    n_treatment = 10000,
    single_arm = FALSE,
    alternative = "less",
    h0 = 0,
    prior_bin = c(1, 1),
    bin_method = "quadrature",
    N_mcmc = NULL
  ))
  # The true probability is approximately 0.7500125; the original answer
  # lies below even the much weaker bound from the posterior moments.
  expect_lt(actual$value$success, 1e-6)
  expect_length(actual$warnings, 1)
  expect_equal(
    actual$warnings[[1]]$counts,
    c(
      flagged_calculations = 1,
      symmetry = 0,
      moments = 1
    )
  )

  lower <- quadrature_beta(1, 3)
  higher <- quadrature_beta(3, 1)
  # Mean difference -1/2 and variance 3/40 imply Pr(D < 0) >= 10/13.
  for (direction in c("less", "greater")) {
    convert <- function(p) if (direction == "less") p else 1 - p
    expect_no_warning(check_beta_quadrature(
      convert(10 / 13 - 5e-7),
      lower,
      higher,
      0,
      direction
    ))
    expect_warning(
      check_beta_quadrature(
        convert(10 / 13 - 2e-6),
        lower,
        higher,
        0,
        direction
      ),
      class = "goldilocks_quadrature_warning"
    )
    # Reversing the arms gives the upper bound 3/13.
    expect_no_warning(check_beta_quadrature(
      convert(0.1),
      higher,
      lower,
      0,
      direction
    ))
    expect_warning(
      check_beta_quadrature(
        convert(0.3),
        higher,
        lower,
        0,
        direction
      ),
      class = "goldilocks_quadrature_warning"
    )
  }
  # A nonzero margin must shift the mean before the bound is calculated.
  uniform <- quadrature_beta(1, 1)
  expect_warning(
    check_beta_quadrature(
      0.5,
      uniform,
      uniform,
      0.5,
      "less"
    ),
    class = "goldilocks_quadrature_warning"
  )
  expect_no_warning(check_beta_quadrature(0.9, uniform, uniform, 0.5, "less"))
  expect_warning(
    check_beta_quadrature(
      0.5,
      uniform,
      uniform,
      -0.5,
      "greater"
    ),
    class = "goldilocks_quadrature_warning"
  )
  expect_no_warning(check_beta_quadrature(
    0.3,
    uniform,
    quadrature_beta(2, 2),
    0,
    "less"
  ))
})

test_that("checks are silent for regular quadrature and other binary methods", {
  args <- list(
    events_control = 15,
    n_control = 50,
    events_treatment = 10,
    n_treatment = 50,
    single_arm = FALSE,
    alternative = "less",
    h0 = 0,
    prior_bin = c(1, 1),
    bin_method = "quadrature",
    N_mcmc = NULL
  )
  expect_no_warning(do.call(bayes_binomial_from_counts, args))
  args$events_control <- args$events_treatment <- 0
  args$n_control <- args$n_treatment <- 10000
  for (method in c("normal", "mc")) {
    args$bin_method <- method
    args$N_mcmc <- 100
    expect_no_warning(do.call(bayes_binomial_from_counts, args))
  }
  args$single_arm <- TRUE
  args$bin_method <- "quadrature"
  expect_no_warning(do.call(bayes_binomial_from_counts, args))
  expect_warning(
    capture_quadrature_warnings(warning("unrelated warning")),
    "unrelated warning"
  )
})

test_that("sim_trials warns once per run and does not alter outcomes or RNG", {
  args <- c(
    quadrature_warning_args(),
    list(Fn = 0, Sn = 1, return_trace = TRUE)
  )
  set.seed(194)
  rng <- .Random.seed
  actual <- collect_quadrature_warnings(do.call(sim_trials, args))
  expect_identical(.Random.seed, rng)
  expect_length(actual$warnings, 1)
  diagnostics <- attr(actual$value, "quadrature_diagnostics")
  expect_gt(diagnostics[["flagged_calculations"]], args$N_trials)
  expect_equal(diagnostics[["absolute_tolerance"]], 1e-6)
  expect_equal(nrow(actual$value$failures), 0)
  expect_equal(actual$warnings[[1]]$counts, diagnostics[1:3])
  expect_identical(
    attr(summarise_sims(actual$value), "quadrature_diagnostics"),
    diagnostics
  )
  summaries <- summarise_sims(list(first = actual$value, second = actual$value))
  expect_identical(
    attr(summaries, "simulation_metadata")[[1]]$quadrature_diagnostics,
    diagnostics
  )
  again <- collect_quadrature_warnings(do.call(sim_trials, args))
  expect_length(again$warnings, 1)
  expect_identical(attr(again$value, "quadrature_diagnostics"), diagnostics)

  # Disable only the checks to recover the pre-warning calculation.
  local_mocked_bindings(check_beta_quadrature = function(...) invisible(NULL))
  expect_no_warning(reference <- do.call(sim_trials, args))
  expect_identical(actual$value$sims, reference$sims)
  expect_identical(actual$value$traces, reference$traces)
  expect_identical(.Random.seed, rng)
  expect_null(attr(reference, "quadrature_diagnostics"))
})

test_that("individual trial and observed-interim calls collect imputation warnings", {
  args <- quadrature_warning_args()
  args[c("N_trials", "seed", "backend")] <- NULL
  args$Fn <- 0
  args$Sn <- 1
  for (traced in c(FALSE, TRUE)) {
    args$return_trace <- traced
    set.seed(441)
    actual <- collect_quadrature_warnings(do.call(survival_adapt, args))
    expect_length(actual$warnings, 1)
    expect_gt(
      attr(actual$value, "quadrature_diagnostics")[["flagged_calculations"]],
      1
    )
  }
  observed <- data.frame(
    id = 1:20,
    treatment = rep(0:1, 10),
    enrollment = 0,
    time = 12,
    event = 0,
    status = "complete"
  )
  actual <- collect_quadrature_warnings(evaluate_interim(
    data = observed,
    data_cut = 12,
    look = 1,
    N_total = 40,
    end_of_study = 12,
    prior_bin = c(1, 10001),
    prior_surv = c(1, 1e10),
    method = "bayes-bin",
    bin_method = "quadrature",
    alternative = "less",
    Fn = 0.05,
    Sn = 1,
    N_impute = 3,
    N_mcmc = NULL,
    seed = 441
  ))
  expect_length(actual$warnings, 1)
  expect_gt(
    attr(actual$value, "quadrature_diagnostics")[["flagged_calculations"]],
    1
  )
})

test_that("saved paths retain diagnostics without warning again on rule application", {
  for (horizons in c(FALSE, TRUE)) {
    args <- quadrature_warning_args()
    if (horizons) {
      args[c("interim_look", "end_of_study")] <- NULL
      args$endpoint_schedule <- data.frame(N = c(20, 30), horizon = c(12, 18))
      args$final_horizon <- 18
    }
    simulate <- if (horizons) sim_horizon_paths else sim_trial_paths
    actual <- collect_quadrature_warnings(do.call(simulate, args))
    expect_length(actual$warnings, 1)
    expect_equal(nrow(actual$value$failures), 0)
    diagnostics <- attr(actual$value, "quadrature_diagnostics")
    counts <- Reduce(
      `+`,
      lapply(actual$value$trials, function(x) {
        attr(x, "quadrature_diagnostics")[1:3]
      })
    )
    expect_equal(diagnostics[1:3], counts)
    if (horizons) {
      expect_no_warning(
        replay <- apply_horizon_rules(
          actual$value,
          args$endpoint_schedule,
          Fn = 0,
          Sn = 1
        )
      )
      direct <- collect_quadrature_warnings(do.call(
        sim_horizon_trials,
        c(args, list(Fn = 0, Sn = 1))
      ))
      expect_length(direct$warnings, 1)
      expect_identical(direct$value$sims, replay$sims)
      expect_gt(
        attr(direct$value, "quadrature_diagnostics")[["flagged_calculations"]],
        1
      )
    } else {
      expect_no_warning(
        replay <- apply_stopping_rules(actual$value, Fn = 0, Sn = 1)
      )
    }
    expect_identical(attr(replay, "quadrature_diagnostics"), diagnostics)
    expect_identical(
      attr(summarise_sims(replay), "quadrature_diagnostics"),
      diagnostics
    )
  }
})

test_that("warnings from final imputations are combined without changing the result", {
  args <- quadrature_warning_args()
  args$interim_look <- NULL
  args$prop_loss <- 0.5
  args$imputed_final <- TRUE
  actual <- collect_quadrature_warnings(do.call(sim_trials, args))
  expect_length(actual$warnings, 1)
  expect_equal(
    attr(actual$value, "quadrature_diagnostics")[["flagged_calculations"]],
    args$N_trials * args$N_impute
  )
  expect_equal(nrow(actual$value$failures), 0)
  local_mocked_bindings(check_beta_quadrature = function(...) invisible(NULL))
  reference <- do.call(sim_trials, args)
  expect_identical(actual$value$sims, reference$sims)
})

test_that("parallel workers return one combined warning to the caller", {
  for (simulate in list(sim_trials, sim_trial_paths, sim_horizon_paths)) {
    args <- quadrature_warning_args()
    if (identical(simulate, sim_trials)) {
      args$Fn <- 0
      args$Sn <- 1
    } else if (identical(simulate, sim_horizon_paths)) {
      args[c("interim_look", "end_of_study")] <- NULL
      args$endpoint_schedule <- data.frame(N = c(20, 30), horizon = c(12, 18))
      args$final_horizon <- 18
    }
    sequential <- collect_quadrature_warnings(do.call(simulate, args))
    for (backend in c("psock", if (.Platform$OS.type != "windows") "fork")) {
      args$backend <- backend
      args$ncores <- 2
      actual <- collect_quadrature_warnings(do.call(simulate, args))
      expect_length(actual$warnings, 1)
      expect_identical(
        attr(actual$value, "quadrature_diagnostics"),
        attr(sequential$value, "quadrature_diagnostics")
      )
      expect_identical(actual$value$failures, sequential$value$failures)
      if (identical(simulate, sim_trials)) {
        expect_identical(actual$value$sims, sequential$value$sims)
      } else {
        expect_identical(actual$value$trials, sequential$value$trials)
      }
    }
  }
})

test_that("a later analysis error does not lose earlier diagnostic warnings", {
  local_mocked_bindings(analyse_final = function(...) {
    stop("final analysis failed")
  })
  args <- c(quadrature_warning_args(), list(Fn = 0, Sn = 1))
  actual <- collect_quadrature_warnings(tryCatch(
    do.call(sim_trials, args),
    error = identity
  ))
  expect_s3_class(actual$value, "goldilocks_all_trials_failed")
  expect_match(conditionMessage(actual$value), "final analysis failed")
  expect_length(actual$warnings, 1)
  expect_gt(
    attr(actual$value, "quadrature_diagnostics")[["flagged_calculations"]],
    1
  )
})
