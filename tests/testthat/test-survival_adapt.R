test_that("survival_adapt-bayes", {
  out <- survival_adapt(
    hazard_treatment = -log(0.85) / 36,
    hazard_control = -log(0.7) / 36,
    cutpoint = 0,
    N_total = 400,
    lambda = 20,
    lambda_time = NULL,
    interim_look = 200,
    end_of_study = 36,
    prior = c(0.1, 0.1),
    block = 2,
    rand_ratio = c(1, 1),
    prop_loss_to_followup = 0.30,
    alternative = "less",
    h0 = 0,
    futility_prob = 0.05,
    expected_success_prob = 0.9,
    prob_ha = 0.975,
    N_impute = 2,
    N_mcmc = 2,
    method = "bayes")

  expect_s3_class(out, "data.frame")
})

test_that("survival_adapt-logrank", {
  out <- survival_adapt(
    hazard_treatment = -log(0.85) / 36,
    hazard_control = -log(0.7) / 36,
    cutpoint = 0,
    N_total = 400,
    lambda = 20,
    lambda_time = NULL,
    interim_look = 200,
    end_of_study = 36,
    prior = c(0.1, 0.1),
    block = 2,
    rand_ratio = c(1, 1),
    prop_loss_to_followup = 0.30,
    alternative = "two.sided",
    h0 = 0,
    futility_prob = 0.05,
    expected_success_prob = 0.9,
    prob_ha = 0.975,
    N_impute = 2,
    N_mcmc = 2,
    method = "logrank")

  expect_s3_class(out, "data.frame")
})

test_that("error-interim-looks", {
  expect_error(
    out <- survival_adapt(
      hazard_treatment = -log(0.85) / 36,
      hazard_control = -log(0.7) / 36,
      cutpoint = 0,
      N_total = 400,
      lambda = 20,
      lambda_time = NULL,
      interim_look = c(200, 300, 400, 500),
      end_of_study = 36,
      prior = c(0.1, 0.1),
      block = 2,
      rand_ratio = c(1, 1),
      prop_loss_to_followup = 0.30,
      alternative = "two.sided",
      h0 = 0,
      futility_prob = 0.05,
      expected_success_prob = 0.9,
      prob_ha = 0.975,
      N_impute = 2,
      N_mcmc = 2,
      method = "logrank")
  )
})

test_that("error-alternative", {
  expect_error(
    out <- survival_adapt(
      hazard_treatment = -log(0.85) / 36,
      hazard_control = -log(0.7) / 36,
      cutpoint = 0,
      N_total = 400,
      lambda = 20,
      lambda_time = NULL,
      interim_look = 200,
      end_of_study = 36,
      prior = c(0.1, 0.1),
      block = 2,
      rand_ratio = c(1, 1),
      prop_loss_to_followup = 0.30,
      alternative = "two-sided",
      h0 = 0,
      futility_prob = 0.05,
      expected_success_prob = 0.9,
      prob_ha = 0.975,
      N_impute = 2,
      N_mcmc = 2,
      method = "logrank")
  )
})

test_that("error-cutpoint", {
  expect_error(
    out <- survival_adapt(
      hazard_treatment = -log(0.85) / 36,
      hazard_control = -log(0.7) / 36,
      cutpoint = 37,
      N_total = 400,
      lambda = 20,
      lambda_time = NULL,
      interim_look = 200,
      end_of_study = 36,
      prior = c(0.1, 0.1),
      block = 2,
      rand_ratio = c(1, 1),
      prop_loss_to_followup = 0.30,
      alternative = "two.sided",
      h0 = 0,
      futility_prob = 0.05,
      expected_success_prob = 0.9,
      prob_ha = 0.975,
      N_impute = 2,
      N_mcmc = 2,
      method = "logrank")
  )
})

test_that("error-alternative-bayes", {
  expect_error(
    out <- survival_adapt(
      hazard_treatment = -log(0.85) / 36,
      hazard_control = -log(0.7) / 36,
      cutpoint = 0,
      N_total = 400,
      lambda = 20,
      lambda_time = NULL,
      interim_look = 200,
      end_of_study = 36,
      prior = c(0.1, 0.1),
      block = 2,
      rand_ratio = c(1, 1),
      prop_loss_to_followup = 0.30,
      alternative = "two.sided",
      h0 = 0,
      futility_prob = 0.05,
      expected_success_prob = 0.9,
      prob_ha = 0.975,
      N_impute = 2,
      N_mcmc = 2,
      method = "bayes")
  )
})

test_that("error-alternative-logrank", {
  expect_error(
    out <- survival_adapt(
      hazard_treatment = -log(0.85) / 36,
      hazard_control = -log(0.7) / 36,
      cutpoint = 0,
      N_total = 400,
      lambda = 20,
      lambda_time = NULL,
      interim_look = 200,
      end_of_study = 36,
      prior = c(0.1, 0.1),
      block = 2,
      rand_ratio = c(1, 1),
      prop_loss_to_followup = 0.30,
      alternative = "less",
      h0 = 0,
      futility_prob = 0.05,
      expected_success_prob = 0.9,
      prob_ha = 0.975,
      N_impute = 2,
      N_mcmc = 2,
      method = "logrank")
  )
})

test_that("error-logrank-single_arm", {
  expect_error(
    out <- survival_adapt(
      hazard_treatment = -log(0.85) / 36,
      hazard_control = NULL,
      cutpoint = 0,
      N_total = 400,
      lambda = 20,
      lambda_time = NULL,
      interim_look = 200,
      end_of_study = 36,
      prior = c(0.1, 0.1),
      block = 2,
      rand_ratio = c(1, 1),
      prop_loss_to_followup = 0.30,
      alternative = "two.sided",
      h0 = 0,
      futility_prob = 0.05,
      expected_success_prob = 0.9,
      prob_ha = 0.975,
      N_impute = 2,
      N_mcmc = 2,
      method = "logrank")
  )
})

test_that("error-prob-thresholds-length_v1", {
  expect_error(
    out <- survival_adapt(
      hazard_treatment = -log(0.85) / 36,
      hazard_control = -log(0.7) / 36,
      cutpoint = 0,
      N_total = 400,
      lambda = 20,
      lambda_time = NULL,
      interim_look = c(100, 200),
      end_of_study = 36,
      prior = c(0.1, 0.1),
      block = 2,
      rand_ratio = c(1, 1),
      prop_loss_to_followup = 0.30,
      alternative = "two.sided",
      h0 = 0,
      futility_prob = 0.05,
      expected_success_prob = c(0.99, 0.9),
      prob_ha = 0.975,
      N_impute = 2,
      N_mcmc = 2,
      method = "logrank")
  )
})

test_that("error-prob-thresholds-length_v2", {
  expect_error(
    out <- survival_adapt(
      hazard_treatment = -log(0.85) / 36,
      hazard_control = -log(0.7) / 36,
      cutpoint = 0,
      N_total = 400,
      lambda = 20,
      lambda_time = NULL,
      interim_look = c(100, 200),
      end_of_study = 36,
      prior = c(0.1, 0.1),
      block = 2,
      rand_ratio = c(1, 1),
      prop_loss_to_followup = 0.30,
      alternative = "two.sided",
      h0 = 0,
      futility_prob = c(0.05, 0.05, 0.05),
      expected_success_prob = c(0.99, 0.99, 0.99),
      prob_ha = 0.975,
      N_impute = 2,
      N_mcmc = 2,
      method = "logrank")
  )
})
