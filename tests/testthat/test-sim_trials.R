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
    bin_N = 100,
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

test_that("sim_trials-zero_cores", {
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
      ncores = ncores,
      seed = 4101
    )
  }

  expect_equal(run_with_cores(1)$sims, run_with_cores(2)$sims)
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
