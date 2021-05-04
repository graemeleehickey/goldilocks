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
    ncores = 1)

  expect_type(out, "list")
  expect_s3_class(out$sims, "data.frame")

  summ_out <- summarise_sims(out$sims)
  expect_s3_class(summ_out, "data.frame")
})

test_that("sim_trials-bayes", {
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
    method = "bayes",
    ncores = 1)

  expect_type(out, "list")
  expect_s3_class(out$sims, "data.frame")

  summ_out <- summarise_sims(out$sims)
  expect_s3_class(summ_out, "data.frame")
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
      method = "bayes",
      ncores = NULL)
  )
})
