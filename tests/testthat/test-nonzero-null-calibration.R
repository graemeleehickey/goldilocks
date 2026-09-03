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

  rejected <- out$sims$post_prob_ha > out$sims$prob_threshold

  expect_mc_rate(
    rejected,
    target = 0.025,
    estimand = "Cox one-sided type-I error"
  )
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

  rejected <- out$sims$post_prob_ha > out$sims$prob_threshold

  expect_mc_rate(
    rejected,
    target = 0.025,
    estimand = "risk-difference one-sided type-I error"
  )
})

test_that("log-rank simulations control one-sided type-I error", {
  out <- sim_trials(
    hazard_treatment = 0.03,
    hazard_control = 0.03,
    N_total = 200,
    lambda = 100,
    end_of_study = 36,
    alternative = "less",
    h0 = 0,
    prob_ha = 0.975,
    N_impute = 1,
    N_mcmc = 1,
    N_trials = 500,
    method = "logrank",
    backend = "sequential",
    seed = 7377
  )

  rejected <- out$sims$post_prob_ha > out$sims$prob_threshold

  expect_mc_rate(
    rejected,
    target = 0.025,
    estimand = "log-rank one-sided type-I error"
  )
})

test_that("Cox Wald confidence intervals have nominal coverage", {
  skip_on_cran()

  set.seed(804)
  n_trials <- 1000L
  true_log_hazard_ratio <- log(0.7)
  critical_value <- stats::qnorm(0.975)
  covered <- replicate(n_trials, {
    data <- sim_comp_data(
      hazard_treatment = 0.03 * 0.7,
      hazard_control = 0.03,
      N_total = 200,
      lambda = 100,
      end_of_study = 36
    )
    fit <- cox_wald_test_checked(data)
    abs(fit$estimate - true_log_hazard_ratio) <= critical_value * fit$std_error
  })

  expect_mc_rate(
    covered,
    target = 0.95,
    estimand = "Cox Wald 95% confidence-interval coverage"
  )
})

test_that("risk-difference Wald confidence intervals have nominal coverage", {
  skip_on_cran()

  set.seed(805)
  n_trials <- 1000L
  control_risk <- 0.2
  treatment_risk <- 0.3
  true_risk_difference <- treatment_risk - control_risk
  end_of_study <- 36
  critical_value <- stats::qnorm(0.975)
  covered <- replicate(n_trials, {
    data <- sim_comp_data(
      hazard_treatment = -log1p(-treatment_risk) / end_of_study,
      hazard_control = -log1p(-control_risk) / end_of_study,
      N_total = 200,
      lambda = 100,
      end_of_study = end_of_study
    )
    fit <- risk_difference_wald_test_checked(
      data = data,
      end_of_study = end_of_study,
      alternative = "greater",
      h0 = 0
    )
    abs(fit$estimate - true_risk_difference) <= critical_value * fit$std_error
  })

  expect_mc_rate(
    covered,
    target = 0.95,
    estimand = "risk-difference Wald 95% confidence-interval coverage"
  )
})
