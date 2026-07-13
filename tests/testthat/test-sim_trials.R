test_that("sim_trials-logrank", {
  hc <- prop_to_haz(c(0.20, 0.30), c(0, 12), 36)
  ht <- prop_to_haz(c(0.05, 0.15), c(0, 12), 36)

  out <- sim_trials(
    hazard_treatment = ht,
    hazard_control = hc,
    cutpoints = c(0, 12),
    N_total = 600,
    lambda = 20,
    lambda_time = 0,
    interim_look = c(400, 500),
    end_of_study = 36,
    prior = c(0.1, 0.1),
    block = 2,
    rand_ratio = c(1, 1),
    prop_loss = 0.30,
    alternative = "two.sided",
    h0 = 0,
    Fn = 0.05,
    Sn = 0.9,
    prob_ha = 0.975,
    N_impute = 5,
    N_mcmc = 2,
    N_trials = 2,
    method = "logrank",
    ncores = 1
  )

  expect_type(out, "list")
  expect_s3_class(out$sims, "data.frame")

  summ_out <- summarise_sims(out$sims)
  expect_s3_class(summ_out, "data.frame")
})

test_that("sim_trials-bayes-surv", {
  hc <- prop_to_haz(c(0.20, 0.30), c(0, 12), 36)
  ht <- prop_to_haz(c(0.05, 0.15), c(0, 12), 36)

  out <- sim_trials(
    hazard_treatment = ht,
    hazard_control = hc,
    cutpoints = c(0, 12),
    N_total = 600,
    lambda = 20,
    lambda_time = 0,
    interim_look = c(400, 500),
    end_of_study = 36,
    prior = c(0.1, 0.1),
    block = 2,
    rand_ratio = c(1, 1),
    prop_loss = 0.30,
    alternative = "less",
    h0 = 0,
    Fn = 0.05,
    Sn = 0.9,
    prob_ha = 0.975,
    N_impute = 5,
    N_mcmc = 2,
    N_trials = 2,
    method = "bayes-surv",
    ncores = 1
  )

  expect_type(out, "list")
  expect_s3_class(out$sims, "data.frame")

  summ_out <- summarise_sims(out$sims)
  expect_s3_class(summ_out, "data.frame")
})

test_that("sim_trials-bayes-bin", {
  hc <- -log(0.7) / 36
  ht <- -log(0.85) / 36

  out <- sim_trials(
    hazard_treatment = ht,
    hazard_control = hc,
    cutpoints = 0,
    N_total = 200,
    lambda = 20,
    lambda_time = 0,
    interim_look = 100,
    end_of_study = 36,
    prior = c(0.1, 0.1),
    bin_prior = c(1, 1),
    bin_method = "normal",
    block = 2,
    rand_ratio = c(1, 1),
    prop_loss = 0.30,
    alternative = "less",
    h0 = 0,
    Fn = 0.05,
    Sn = 0.9,
    prob_ha = 0.975,
    N_impute = 2,
    N_mcmc = 2,
    N_trials = 1,
    method = "bayes-bin",
    ncores = 1,
    seed = 9241
  )

  expect_type(out, "list")
  expect_s3_class(out$sims, "data.frame")
})

test_that("sim_trials rejects invalid ncores", {
  hc <- prop_to_haz(c(0.20, 0.30), c(0, 12), 36)
  ht <- prop_to_haz(c(0.05, 0.15), c(0, 12), 36)

  expect_error(
    sim_trials(
      hazard_treatment = ht,
      hazard_control = hc,
      cutpoints = c(0, 12),
      N_total = 600,
      lambda = 20,
      lambda_time = 0,
      interim_look = c(400, 500),
      end_of_study = 36,
      prior = c(0.1, 0.1),
      block = 2,
      rand_ratio = c(1, 1),
      prop_loss = 0.30,
      alternative = "less",
      h0 = 0,
      Fn = 0.05,
      Sn = 0.9,
      prob_ha = 0.975,
      N_impute = 5,
      N_mcmc = 2,
      N_trials = 2,
      method = "bayes-surv",
      ncores = NULL
    )
  )
})

test_that("sim_trials defaults to serial execution", {
  expect_identical(formals(sim_trials)$ncores, 1L)
  expect_identical(formals(sim_trials)$return_trace, FALSE)
  expect_equal(
    as.character(formals(sim_trials)$backend)[-1],
    c("auto", "fork", "psock", "sequential")
  )
})

test_that("sim_trials optionally retains traces without changing summaries", {
  args <- list(
    hazard_treatment = -log(0.85) / 36,
    hazard_control = -log(0.7) / 36,
    cutpoints = 0,
    N_total = 80,
    lambda = 20,
    lambda_time = 0,
    interim_look = c(40, 60),
    end_of_study = 36,
    prior = c(0.1, 0.1),
    block = 2,
    rand_ratio = c(1, 1),
    prop_loss = 0.05,
    alternative = "less",
    h0 = 0,
    Fn = c(0.05, 0.05),
    Sn = c(0.95, 0.9),
    prob_ha = 0.95,
    N_impute = 2,
    N_mcmc = 2,
    N_trials = 2,
    method = "bayes-surv",
    ncores = 1,
    seed = 4102
  )

  compact <- do.call(sim_trials, args)
  traced <- do.call(sim_trials, c(args, list(return_trace = TRUE)))

  expect_equal(traced$sims, compact$sims)
  expect_named(traced, c("sims", "traces", "call"))
  expect_s3_class(traced$traces, "data.frame")
  expect_true(all(c(
    "trial", "look", "ppp_stop_now", "ppp_success_at_max", "decision"
  ) %in% names(traced$traces)))
  expect_setequal(unique(traced$traces$trial), 1:2)
})

test_that("sim_trials validates return_trace", {
  expect_error(
    sim_trials(
      hazard_treatment = -log(0.85) / 36,
      hazard_control = -log(0.7) / 36,
      N_total = 40,
      interim_look = 20,
      end_of_study = 36,
      return_trace = NA
    ),
    "return_trace"
  )
})

test_that("sim_trials uses survival_adapt decision-rule defaults", {
  expect_identical(
    formals(sim_trials)$alternative,
    formals(survival_adapt)$alternative
  )
  expect_identical(formals(sim_trials)$Fn, formals(survival_adapt)$Fn)
})

test_that("sim_trials validates N_trials", {
  hc <- -log(0.7) / 36
  ht <- -log(0.85) / 36

  expect_error(
    sim_trials(
      hazard_treatment = ht,
      hazard_control = hc,
      cutpoints = 0,
      N_total = 50,
      lambda = 20,
      lambda_time = 0,
      interim_look = NULL,
      end_of_study = 36,
      N_trials = 1.5,
      method = "logrank",
      ncores = 1
    ),
    "N_trials"
  )
})

test_that("sim_trials is reproducible with seed and ncores = 1", {
  hc <- -log(0.7) / 36
  ht <- -log(0.85) / 36

  run_once <- function() {
    sim_trials(
      hazard_treatment = ht,
      hazard_control = hc,
      cutpoints = 0,
      N_total = 200,
      lambda = 20,
      lambda_time = 0,
      interim_look = 100,
      end_of_study = 36,
      prior = c(0.1, 0.1),
      block = 2,
      rand_ratio = c(1, 1),
      prop_loss = 0.30,
      alternative = "two.sided",
      h0 = 0,
      Fn = 0.05,
      Sn = 0.9,
      prob_ha = 0.975,
      N_impute = 2,
      N_mcmc = 2,
      N_trials = 2,
      method = "logrank",
      ncores = 1,
      seed = 4101
    )
  }

  expect_equal(run_once()$sims, run_once()$sims)
})

test_that("sim_trials uses reproducible per-trial streams in parallel", {
  skip_on_os("windows")

  hc <- -log(0.7) / 36
  ht <- -log(0.85) / 36

  run_with_cores <- function(ncores) {
    sim_trials(
      hazard_treatment = ht,
      hazard_control = hc,
      cutpoints = 0,
      N_total = 200,
      lambda = 20,
      lambda_time = 0,
      interim_look = 100,
      end_of_study = 36,
      prior = c(0.1, 0.1),
      block = 2,
      rand_ratio = c(1, 1),
      prop_loss = 0.30,
      alternative = "two.sided",
      h0 = 0,
      Fn = 0.05,
      Sn = 0.9,
      prob_ha = 0.975,
      N_impute = 2,
      N_mcmc = 2,
      N_trials = 2,
      method = "logrank",
      return_trace = TRUE,
      ncores = ncores,
      seed = 4101
    )
  }

  sequential <- run_with_cores(1)
  parallel <- run_with_cores(2)
  expect_equal(sequential$sims, parallel$sims)
  expect_equal(sequential$traces, parallel$traces)
})

test_that("sim_trials produces identical seeded PSOCK results", {
  hc <- -log(0.7) / 36
  ht <- -log(0.85) / 36

  run_with_backend <- function(backend, ncores) {
    sim_trials(
      hazard_treatment = ht,
      hazard_control = hc,
      cutpoints = 0,
      N_total = 100,
      lambda = 20,
      lambda_time = 0,
      interim_look = 50,
      end_of_study = 36,
      prior = c(0.1, 0.1),
      block = 2,
      rand_ratio = c(1, 1),
      prop_loss = 0.30,
      alternative = "two.sided",
      h0 = 0,
      Fn = 0.05,
      Sn = 0.9,
      prob_ha = 0.975,
      N_impute = 2,
      N_mcmc = 2,
      N_trials = 2,
      method = "logrank",
      return_trace = TRUE,
      ncores = ncores,
      backend = backend,
      seed = 4101
    )
  }

  sequential <- run_with_backend("sequential", 1)
  psock <- run_with_backend("psock", 2)
  expect_equal(sequential$sims, psock$sims)
  expect_equal(sequential$traces, psock$traces)
})

test_that("sim_trials auto backend selects the platform default", {
  expected <- if (.Platform$OS.type == "windows") "psock" else "fork"
  expect_identical(resolve_sim_backend("auto", 2L), expected)
  expect_identical(resolve_sim_backend("auto", 1L), "sequential")
})

test_that("sim_trials preserves the caller RNG state when seeded", {
  set.seed(9012)
  before_kind <- RNGkind()
  before <- .Random.seed

  sim_trials(
    hazard_treatment = -log(0.85) / 36,
    hazard_control = -log(0.7) / 36,
    cutpoints = 0,
    N_total = 50,
    lambda = 20,
    lambda_time = 0,
    interim_look = NULL,
    end_of_study = 36,
    N_trials = 1,
    method = "logrank",
    ncores = 1,
    seed = 4101
  )

  expect_identical(RNGkind(), before_kind)
  expect_identical(.Random.seed, before)
})

test_that("sim_trials validates seed", {
  hc <- -log(0.7) / 36
  ht <- -log(0.85) / 36

  expect_error(
    sim_trials(
      hazard_treatment = ht,
      hazard_control = hc,
      cutpoints = 0,
      N_total = 200,
      lambda = 20,
      lambda_time = 0,
      interim_look = 100,
      end_of_study = 36,
      prior = c(0.1, 0.1),
      block = 2,
      rand_ratio = c(1, 1),
      prop_loss = 0.30,
      alternative = "two.sided",
      h0 = 0,
      Fn = 0.05,
      Sn = 0.9,
      prob_ha = 0.975,
      N_impute = 2,
      N_mcmc = 2,
      N_trials = 2,
      method = "logrank",
      ncores = 1,
      seed = c(1, 2)
    ),
    "single integer"
  )
})
