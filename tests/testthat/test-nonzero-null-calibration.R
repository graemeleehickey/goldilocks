test_that("Cox nonzero-null simulations control one-sided type-I error", {
  null_margin <- 1.5
  out <- sim_trials(
    hazard_treatment = 0.03 * null_margin,
    hazard_control = 0.03,
    N_total = 200,
    lambda = 100,
    end_of_study = 36,
    alternative = "less",
    h0 = log(null_margin),
    prob_ha = 0.975,
    N_impute = 1,
    N_mcmc = 1,
    N_trials = 500,
    method = "cox",
    backend = "sequential",
    seed = 7373
  )

  rejection_rate <- mean(
    out$sims$post_prob_ha > out$sims$prob_threshold
  )
  expect_gt(rejection_rate, 0)
  expect_lte(rejection_rate, 0.06)
})

test_that("risk-difference nonzero-null simulations control type-I error", {
  control_risk <- 0.2
  null_difference <- 0.1
  out <- sim_trials(
    hazard_treatment = -log1p(-(control_risk + null_difference)) / 36,
    hazard_control = -log1p(-control_risk) / 36,
    N_total = 200,
    lambda = 100,
    end_of_study = 36,
    alternative = "greater",
    h0 = null_difference,
    prob_ha = 0.975,
    N_impute = 1,
    N_mcmc = 1,
    N_trials = 500,
    method = "riskdiff",
    backend = "sequential",
    seed = 7375
  )

  rejection_rate <- mean(
    out$sims$post_prob_ha > out$sims$prob_threshold
  )
  expect_gt(rejection_rate, 0)
  expect_lte(rejection_rate, 0.06)
})
