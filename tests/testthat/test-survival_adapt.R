test_that("survival_adapt-bayes", {
  out <- survival_adapt(
    hazard_treatment = -log(0.85) / 36,
    hazard_control = -log(0.7) / 36,
    cutpoints = 0,
    N_total = 400,
    lambda = 20,
    lambda_time = 0,
    interim_look = 200,
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
    N_impute = 2,
    N_mcmc = 2,
    method = "bayes")

  expect_s3_class(out, "data.frame")
})

test_that("survival_adapt-logrank", {
  out <- survival_adapt(
    hazard_treatment = -log(0.85) / 36,
    hazard_control = -log(0.7) / 36,
    cutpoints = 0,
    N_total = 400,
    lambda = 20,
    lambda_time = 0,
    interim_look = 200,
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
    method = "logrank")

  expect_s3_class(out, "data.frame")
})

test_that("survival_adapt-cox", {
  out <- survival_adapt(
    hazard_treatment = -log(0.85) / 36,
    hazard_control = -log(0.7) / 36,
    cutpoints = 0,
    N_total = 400,
    lambda = 20,
    lambda_time = 0,
    interim_look = 200,
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
    method = "cox")

  expect_s3_class(out, "data.frame")
  expect_true(!is.na(out$est_final))
})

test_that("survival_adapt-chisq", {
  set.seed(1)
  out <- survival_adapt(
    hazard_treatment = -log(0.85) / 36,
    hazard_control = -log(0.7) / 36,
    cutpoints = 0,
    N_total = 400,
    lambda = 20,
    lambda_time = 0,
    interim_look = 200,
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
    method = "chisq")

  expect_s3_class(out, "data.frame")
  expect_true(!is.na(out$est_final))
})

test_that("survival_adapt-chisq excludes LTFU when imputed_final = FALSE", {
  set.seed(3927)
  out <- survival_adapt(
    hazard_treatment = -log(0.85) / 36,
    hazard_control = -log(0.7) / 36,
    cutpoints = 0,
    N_total = 400,
    lambda = 20,
    lambda_time = 0,
    interim_look = 200,
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
    method = "chisq",
    imputed_final = FALSE)

  expect_s3_class(out, "data.frame")
  expect_true(out$post_prob_ha >= 0 && out$post_prob_ha <= 1)
})

test_that("survival_adapt-cox with imputed_final", {
  set.seed(8263)
  out <- survival_adapt(
    hazard_treatment = -log(0.85) / 36,
    hazard_control = -log(0.7) / 36,
    cutpoints = 0,
    N_total = 300,
    lambda = 20,
    lambda_time = 0,
    interim_look = 150,
    end_of_study = 36,
    prior = c(0.1, 0.1),
    block = 2,
    rand_ratio = c(1, 1),
    prop_loss = 0,
    alternative = "two.sided",
    h0 = 0,
    Fn = 0.05,
    Sn = 0.9,
    prob_ha = 0.975,
    N_impute = 2,
    N_mcmc = 2,
    method = "cox",
    imputed_final = TRUE)

  expect_s3_class(out, "data.frame")
  expect_true(!is.na(out$post_prob_ha))
})

test_that("survival_adapt-complex", {

  skip_on_cran()

  hc <- prop_to_haz(c(0.20, 0.25, 0.30), c(0, 3, 6), 24)
  ht <- prop_to_haz(c(0.075, 0.125, 0.175), c(0, 3, 6), 24)

  set.seed(12345)
  out <- survival_adapt(
    hazard_treatment = ht,
    hazard_control = hc,
    cutpoints = c(0, 3, 6),
    N_total = 500,
    lambda = c(6, 8),
    lambda_time = c(0, 6),
    interim_look = c(200, 350),
    end_of_study = 24,
    prior = c(0.1, 0.1),
    block = 3,
    rand_ratio = c(2, 1),
    prop_loss = 0,
    alternative = "less",
    h0 = 0,
    Fn = 0.05,
    Sn = 0.9,
    prob_ha = 0.95,
    N_impute = 20,
    N_mcmc = 20,
    method = "bayes",
    imputed_final = TRUE)

  expect_s3_class(out, "data.frame")
})

test_that("error-interim-looks", {
  expect_error(
    out <- survival_adapt(
      hazard_treatment = -log(0.85) / 36,
      hazard_control = -log(0.7) / 36,
      cutpoints = 0,
      N_total = 400,
      lambda = 20,
      lambda_time = 0,
      interim_look = c(200, 300, 400, 500),
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
      method = "logrank")
  )
})

test_that("error-interim-look-below-block-size", {
  # A two-arm interim look smaller than the block size can enrol a single arm
  # only, leaving the interim posterior undefined for the missing arm. This
  # should be caught as an input error before any simulation runs.
  expect_error(
    out <- survival_adapt(
      hazard_treatment = -log(0.85) / 36,
      hazard_control = -log(0.7) / 36,
      cutpoints = 0,
      N_total = 100,
      lambda = 20,
      lambda_time = 0,
      interim_look = 3,
      end_of_study = 36,
      prior = c(0.1, 0.1),
      block = 4,
      rand_ratio = c(1, 1),
      prop_loss = 0.30,
      alternative = "less",
      h0 = 0,
      Fn = 0.05,
      Sn = 0.9,
      prob_ha = 0.975,
      N_impute = 2,
      N_mcmc = 2,
      method = "bayes"),
    "must be at least the block size"
  )
})

test_that("error-alternative", {
  expect_error(
    out <- survival_adapt(
      hazard_treatment = -log(0.85) / 36,
      hazard_control = -log(0.7) / 36,
      cutpoints = 0,
      N_total = 400,
      lambda = 20,
      lambda_time = 0,
      interim_look = 200,
      end_of_study = 36,
      prior = c(0.1, 0.1),
      block = 2,
      rand_ratio = c(1, 1),
      prop_loss = 0.30,
      alternative = "two-sided",
      h0 = 0,
      Fn = 0.05,
      Sn = 0.9,
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
      cutpoints = 37,
      N_total = 400,
      lambda = 20,
      lambda_time = 0,
      interim_look = 200,
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
      method = "logrank")
  )
})

test_that("error-alternative-bayes", {
  expect_error(
    out <- survival_adapt(
      hazard_treatment = -log(0.85) / 36,
      hazard_control = -log(0.7) / 36,
      cutpoints = 0,
      N_total = 400,
      lambda = 20,
      lambda_time = 0,
      interim_look = 200,
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
      method = "bayes")
  )
})

test_that("survival_adapt-logrank-one-sided", {
  out <- survival_adapt(
    hazard_treatment = -log(0.85) / 36,
    hazard_control = -log(0.7) / 36,
    cutpoints = 0,
    N_total = 400,
    lambda = 20,
    lambda_time = 0,
    interim_look = 200,
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
    N_impute = 2,
    N_mcmc = 2,
    method = "logrank")

  expect_s3_class(out, "data.frame")
  expect_true(out$post_prob_ha >= 0 && out$post_prob_ha <= 1)
})

test_that("survival_adapt-cox-one-sided-less", {
  out <- survival_adapt(
    hazard_treatment = -log(0.85) / 36,
    hazard_control = -log(0.7) / 36,
    cutpoints = 0,
    N_total = 400,
    lambda = 20,
    lambda_time = 0,
    interim_look = 200,
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
    N_impute = 2,
    N_mcmc = 2,
    method = "cox")

  expect_s3_class(out, "data.frame")
  expect_true(out$post_prob_ha >= 0 && out$post_prob_ha <= 1)
  expect_true(!is.na(out$est_final))
})

test_that("survival_adapt-cox-one-sided-greater", {
  set.seed(1)
  out <- survival_adapt(
    hazard_treatment = -log(0.85) / 36,
    hazard_control = -log(0.7) / 36,
    cutpoints = 0,
    N_total = 400,
    lambda = 20,
    lambda_time = 0,
    interim_look = 200,
    end_of_study = 36,
    prior = c(0.1, 0.1),
    block = 2,
    rand_ratio = c(1, 1),
    prop_loss = 0.30,
    alternative = "greater",
    h0 = 0,
    Fn = 0.05,
    Sn = 0.9,
    prob_ha = 0.975,
    N_impute = 2,
    N_mcmc = 2,
    method = "cox")

  expect_s3_class(out, "data.frame")
  expect_true(out$post_prob_ha >= 0 && out$post_prob_ha <= 1)
})

test_that("error-alternative-chisq", {
  expect_error(
    survival_adapt(
      hazard_treatment = -log(0.85) / 36,
      hazard_control = -log(0.7) / 36,
      cutpoints = 0,
      N_total = 400,
      lambda = 20,
      lambda_time = 0,
      interim_look = 200,
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
      N_impute = 2,
      N_mcmc = 2,
      method = "chisq"),
    "chi-square"
  )
})

test_that("error-logrank-single_arm", {
  expect_error(
    out <- survival_adapt(
      hazard_treatment = -log(0.85) / 36,
      hazard_control = NULL,
      cutpoints = 0,
      N_total = 400,
      lambda = 20,
      lambda_time = 0,
      interim_look = 200,
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
      method = "logrank")
  )
})

test_that("error-prob-thresholds-length_v1", {
  expect_error(
    out <- survival_adapt(
      hazard_treatment = -log(0.85) / 36,
      hazard_control = -log(0.7) / 36,
      cutpoints = 0,
      N_total = 400,
      lambda = 20,
      lambda_time = 0,
      interim_look = c(100, 200),
      end_of_study = 36,
      prior = c(0.1, 0.1),
      block = 2,
      rand_ratio = c(1, 1),
      prop_loss = 0.30,
      alternative = "two.sided",
      h0 = 0,
      Fn = 0.05,
      Sn = c(0.99, 0.9),
      prob_ha = 0.975,
      N_impute = 2,
      N_mcmc = 2,
      method = "logrank")
  )
})

test_that("survival_adapt works with no interim looks", {
  set.seed(4817)
  out <- survival_adapt(
    hazard_treatment = -log(0.85) / 36,
    hazard_control = -log(0.7) / 36,
    cutpoints = 0,
    N_total = 200,
    lambda = 20,
    lambda_time = 0,
    interim_look = NULL,
    end_of_study = 36,
    prior = c(0.1, 0.1),
    block = 2,
    rand_ratio = c(1, 1),
    prop_loss = 0,
    alternative = "two.sided",
    h0 = 0,
    Fn = NULL,
    Sn = NULL,
    prob_ha = 0.975,
    N_impute = 2,
    N_mcmc = 2,
    method = "logrank")

  expect_s3_class(out, "data.frame")
  expect_equal(out$N_enrolled, 200)
  expect_equal(out$stop_futility, 0)
  expect_equal(out$stop_expected_success, 0)
  expect_true(is.na(out$ppp_success))
})

test_that("survival_adapt keeps futility disabled when Fn = 0", {
  set.seed(3167)
  out <- survival_adapt(
    hazard_treatment = -log(0.85) / 36,
    hazard_control = -log(0.7) / 36,
    cutpoints = 0,
    N_total = 200,
    lambda = 20,
    lambda_time = 0,
    interim_look = c(100, 150),
    end_of_study = 36,
    prior = c(0.1, 0.1),
    block = 2,
    rand_ratio = c(1, 1),
    prop_loss = 0.30,
    alternative = "two.sided",
    h0 = 0,
    Fn = c(0, 0),
    Sn = c(1, 1),
    prob_ha = 0.975,
    N_impute = 2,
    N_mcmc = 2,
    method = "logrank")

  expect_s3_class(out, "data.frame")
  expect_equal(out$stop_futility, 0)
})

test_that("error-prob-thresholds-length_v2", {
  expect_error(
    out <- survival_adapt(
      hazard_treatment = -log(0.85) / 36,
      hazard_control = -log(0.7) / 36,
      cutpoints = 0,
      N_total = 400,
      lambda = 20,
      lambda_time = 0,
      interim_look = c(100, 200),
      end_of_study = 36,
      prior = c(0.1, 0.1),
      block = 2,
      rand_ratio = c(1, 1),
      prop_loss = 0.30,
      alternative = "two.sided",
      h0 = 0,
      Fn = c(0.05, 0.05, 0.05),
      Sn = c(0.99, 0.99, 0.99),
      prob_ha = 0.975,
      N_impute = 2,
      N_mcmc = 2,
      method = "logrank")
  )
})
