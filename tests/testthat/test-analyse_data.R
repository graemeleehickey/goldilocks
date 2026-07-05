test_that("analyse_data works with method = 'logrank'", {
  set.seed(6184)
  data <- data.frame(
    time = rexp(100, 0.02),
    event = sample(0:1, 100, replace = TRUE),
    treatment = rep(0:1, each = 50)
  )
  res <- analyse_data(
    data = data,
    cutpoints = 0,
    end_of_study = 36,
    prior = c(0.1, 0.1),
    N_mcmc = 10,
    single_arm = FALSE,
    method = "logrank",
    alternative = "greater",
    h0 = 0
  )
  expect_type(res, "list")
  expect_named(res, c("success", "effect"))
  expect_true(res$success >= 0 && res$success <= 1)
  expect_true(is.na(res$effect))
})

test_that("analyse_data works with method = 'cox'", {
  set.seed(3729)
  data <- data.frame(
    time = rexp(100, 0.02),
    event = sample(0:1, 100, replace = TRUE),
    treatment = rep(0:1, each = 50)
  )
  res <- analyse_data(
    data = data,
    cutpoints = 0,
    end_of_study = 36,
    prior = c(0.1, 0.1),
    N_mcmc = 10,
    single_arm = FALSE,
    method = "cox",
    alternative = "greater",
    h0 = 0
  )
  expect_type(res, "list")
  expect_named(res, c("success", "effect"))
  expect_true(res$success >= 0 && res$success <= 1)
  expect_type(res$effect, "double")
})

test_that("cox_wald_test matches coxph treatment estimate and standard error", {
  set.seed(4291)
  data <- data.frame(
    time = c(rexp(50, rate = 0.06), rexp(50, rate = 0.04)),
    event = sample(0:1, 100, replace = TRUE),
    treatment = rep(0:1, each = 50)
  )

  fast_fit <- cox_wald_test(data)
  survival_fit <- coxph(Surv(time, event) ~ treatment, data = data)

  expect_equal(fast_fit$estimate, unname(survival_fit$coefficients[1]))
  expect_equal(fast_fit$std_error, sqrt(unname(survival_fit$var[1, 1])))
})

test_that("analyse_data works with method = 'bayes' (two-arm)", {
  set.seed(8415)
  data <- data.frame(
    time = rexp(100, 0.02),
    event = sample(0:1, 100, replace = TRUE),
    treatment = rep(0:1, each = 50)
  )
  res <- analyse_data(
    data = data,
    cutpoints = 0,
    end_of_study = 36,
    prior = c(0.1, 0.1),
    N_mcmc = 100,
    single_arm = FALSE,
    method = "bayes",
    alternative = "greater",
    h0 = 0
  )
  expect_type(res, "list")
  expect_named(res, c("success", "effect"))
  expect_true(res$success >= 0 && res$success <= 1)
  expect_length(res$effect, 100)
})

test_that("analyse_data works with method = 'bayes' and alternative = 'less'", {
  set.seed(5091)
  data <- data.frame(
    time = rexp(100, 0.02),
    event = sample(0:1, 100, replace = TRUE),
    treatment = rep(0:1, each = 50)
  )
  res_greater <- analyse_data(
    data = data,
    cutpoints = 0,
    end_of_study = 36,
    prior = c(0.1, 0.1),
    N_mcmc = 500,
    single_arm = FALSE,
    method = "bayes",
    alternative = "greater",
    h0 = 0
  )
  res_less <- analyse_data(
    data = data,
    cutpoints = 0,
    end_of_study = 36,
    prior = c(0.1, 0.1),
    N_mcmc = 500,
    single_arm = FALSE,
    method = "bayes",
    alternative = "less",
    h0 = 0
  )
  # P(effect > 0) + P(effect < 0) should be approximately 1
  # (exact only if no draws == 0)
  expect_equal(res_greater$success + res_less$success, 1, tolerance = 0.05)
})

test_that("analyse_data works with method = 'chisq'", {
  set.seed(2647)
  event <- sample(0:1, 100, replace = TRUE, prob = c(0.3, 0.7))
  data <- data.frame(
    time = ifelse(event == 1, rexp(100, 0.02), 36),
    event = event,
    treatment = rep(0:1, each = 50)
  )
  res <- analyse_data(
    data = data,
    cutpoints = 0,
    end_of_study = 36,
    prior = c(0.1, 0.1),
    N_mcmc = 10,
    single_arm = FALSE,
    method = "chisq",
    alternative = "greater",
    h0 = 0
  )
  expect_type(res, "list")
  expect_named(res, c("success", "effect"))
  expect_true(res$success >= 0 && res$success <= 1)
  expect_type(res$effect, "double")
})

test_that("analyse_data method = 'chisq' rejects incomplete censored outcomes", {
  data <- data.frame(
    time = c(5, 36, 7, 36),
    event = c(0, 0, 1, 1),
    treatment = c(0, 0, 1, 1)
  )

  expect_error(
    analyse_data(
      data = data,
      cutpoints = 0,
      end_of_study = 36,
      prior = c(0.1, 0.1),
      N_mcmc = 10,
      single_arm = FALSE,
      method = "chisq",
      alternative = "two.sided",
      h0 = 0
    ),
    "requires all censored subjects"
  )
})

test_that("analyse_data works with method = 'bayes-bin' and Monte Carlo", {
  set.seed(4912)
  data <- data.frame(
    time = rep(36, 100),
    event = c(rep(1, 10), rep(0, 40), rep(1, 30), rep(0, 20)),
    treatment = rep(0:1, each = 50)
  )

  res <- analyse_data(
    data = data,
    cutpoints = 0,
    end_of_study = 36,
    prior = c(0.1, 0.1),
    N_mcmc = 10,
    single_arm = FALSE,
    method = "bayes-bin",
    alternative = "greater",
    h0 = 0,
    bin_prior = c(1, 1),
    bin_method = "mc",
    bin_N = 5000
  )

  expect_type(res, "list")
  expect_named(res, c("success", "effect"))
  expect_true(res$success >= 0 && res$success <= 1)
  expect_length(res$effect, 5000)
  expect_true(res$success > 0.9)
})

test_that("analyse_data method = 'bayes-bin' quadrature tails sum to 1", {
  data <- data.frame(
    time = rep(36, 40),
    event = c(rep(1, 8), rep(0, 12), rep(1, 13), rep(0, 7)),
    treatment = rep(0:1, each = 20)
  )
  args <- list(
    data = data,
    cutpoints = 0,
    end_of_study = 36,
    prior = c(0.1, 0.1),
    N_mcmc = 10,
    single_arm = FALSE,
    method = "bayes-bin",
    h0 = 0,
    bin_prior = c(1, 1),
    bin_method = "quadrature",
    bin_N = 1000
  )

  res_less <- do.call(analyse_data, c(args, alternative = "less"))
  res_greater <- do.call(analyse_data, c(args, alternative = "greater"))

  expect_equal(res_less$success + res_greater$success, 1, tolerance = 1e-8)
  expect_type(res_less$effect, "double")
})

test_that("analyse_data method = 'bayes-bin' engines agree on stable data", {
  data <- data.frame(
    time = rep(36, 400),
    event = c(rep(1, 80), rep(0, 120), rep(1, 120), rep(0, 80)),
    treatment = rep(0:1, each = 200)
  )
  args <- list(
    data = data,
    cutpoints = 0,
    end_of_study = 36,
    prior = c(0.1, 0.1),
    N_mcmc = 10,
    single_arm = FALSE,
    method = "bayes-bin",
    alternative = "greater",
    h0 = 0,
    bin_prior = c(1, 1)
  )

  set.seed(3319)
  res_mc <- do.call(analyse_data, c(args, bin_method = "mc", bin_N = 50000))
  res_normal <- do.call(
    analyse_data,
    c(args, bin_method = "normal", bin_N = 1000)
  )
  res_quadrature <- do.call(
    analyse_data,
    c(args, bin_method = "quadrature", bin_N = 1000)
  )

  expect_equal(res_mc$success, res_quadrature$success, tolerance = 0.01)
  expect_equal(res_normal$success, res_quadrature$success, tolerance = 0.02)
})

test_that("analyse_data method = 'bayes-bin' works for single-arm data", {
  data <- data.frame(
    time = rep(36, 50),
    event = c(rep(1, 15), rep(0, 35)),
    treatment = rep(1, 50)
  )

  res <- analyse_data(
    data = data,
    cutpoints = 0,
    end_of_study = 36,
    prior = c(0.1, 0.1),
    N_mcmc = 10,
    single_arm = TRUE,
    method = "bayes-bin",
    alternative = "less",
    h0 = 0.4,
    bin_prior = c(1, 1),
    bin_method = "quadrature",
    bin_N = 1000
  )

  expect_true(res$success > 0.9)
  expect_type(res$effect, "double")
})

test_that("analyse_data method = 'bayes-bin' rejects incomplete censored outcomes", {
  data <- data.frame(
    time = c(36, 5, 36, 36),
    event = c(1, 0, 1, 0),
    treatment = c(0, 0, 1, 1)
  )

  expect_error(
    analyse_data(
      data = data,
      cutpoints = 0,
      end_of_study = 36,
      prior = c(0.1, 0.1),
      N_mcmc = 10,
      single_arm = FALSE,
      method = "bayes-bin",
      alternative = "greater",
      h0 = 0,
      bin_prior = c(1, 1),
      bin_method = "mc",
      bin_N = 1000
    ),
    "requires all censored subjects"
  )
})

test_that("analyse_data works with method = 'bayes' (single-arm)", {
  set.seed(9473)
  data <- data.frame(
    time = rexp(50, 0.02),
    event = sample(0:1, 50, replace = TRUE),
    treatment = rep(1, 50)
  )
  res <- analyse_data(
    data = data,
    cutpoints = 0,
    end_of_study = 36,
    prior = c(0.1, 0.1),
    N_mcmc = 100,
    single_arm = TRUE,
    method = "bayes",
    alternative = "greater",
    h0 = 0
  )
  expect_type(res, "list")
  expect_true(res$success >= 0 && res$success <= 1)
  expect_length(res$effect, 100)
})

test_that("analyse_data works with piecewise cutpoints", {
  set.seed(4318)
  data <- data.frame(
    time = rexp(100, 0.02),
    event = sample(0:1, 100, replace = TRUE),
    treatment = rep(0:1, each = 50)
  )
  res <- analyse_data(
    data = data,
    cutpoints = c(0, 12),
    end_of_study = 36,
    prior = c(0.1, 0.1),
    N_mcmc = 100,
    single_arm = FALSE,
    method = "bayes",
    alternative = "greater",
    h0 = 0
  )
  expect_type(res, "list")
  expect_true(res$success >= 0 && res$success <= 1)
})

test_that("analyse_data one-sided cox: less + greater sum to ~1", {
  set.seed(5813)
  data <- data.frame(
    time = rexp(100, 0.02),
    event = sample(0:1, 100, replace = TRUE),
    treatment = rep(0:1, each = 50)
  )
  args <- list(
    data = data,
    cutpoints = 0,
    end_of_study = 36,
    prior = c(0.1, 0.1),
    N_mcmc = 10,
    single_arm = FALSE,
    method = "cox",
    h0 = 0
  )
  res_less <- do.call(analyse_data, c(args, alternative = "less"))
  res_greater <- do.call(analyse_data, c(args, alternative = "greater"))
  res_two <- do.call(analyse_data, c(args, alternative = "two.sided"))
  # One-sided p-values should sum to 1
  # success = 1 - p, so (1 - success_less) + (1 - success_greater) = 1
  expect_equal(
    (1 - res_less$success) + (1 - res_greater$success),
    1,
    tolerance = 1e-10
  )
  # Two-sided success should be <= max of one-sided
  expect_true(res_two$success <= max(res_less$success, res_greater$success))
})

test_that("analyse_data one-sided logrank: less + greater sum to ~1", {
  set.seed(7249)
  data <- data.frame(
    time = rexp(100, 0.02),
    event = sample(0:1, 100, replace = TRUE),
    treatment = rep(0:1, each = 50)
  )
  args <- list(
    data = data,
    cutpoints = 0,
    end_of_study = 36,
    prior = c(0.1, 0.1),
    N_mcmc = 10,
    single_arm = FALSE,
    method = "logrank",
    h0 = 0
  )
  res_less <- do.call(analyse_data, c(args, alternative = "less"))
  res_greater <- do.call(analyse_data, c(args, alternative = "greater"))
  expect_equal(
    (1 - res_less$success) + (1 - res_greater$success),
    1,
    tolerance = 1e-10
  )
})

test_that("analyse_data one-sided cox detects beneficial treatment", {
  set.seed(8164)
  # Treatment has much lower hazard (beneficial)
  data <- data.frame(
    time = c(rexp(50, rate = 0.10), rexp(50, rate = 0.02)),
    event = rep(1L, 100),
    treatment = rep(0:1, each = 50)
  )
  res <- analyse_data(
    data = data,
    cutpoints = 0,
    end_of_study = 36,
    prior = c(0.1, 0.1),
    N_mcmc = 10,
    single_arm = FALSE,
    method = "cox",
    alternative = "less",
    h0 = 0
  )
  # "less" = treatment beneficial; should have high success
  expect_true(res$success > 0.95)
})

test_that("analyse_data one-sided cox uses h0 as the null log hazard ratio", {
  set.seed(123)
  data <- data.frame(
    time = c(rexp(200, rate = 0.10), rexp(200, rate = 0.12)),
    event = rep(1L, 400),
    treatment = rep(0:1, each = 200)
  )

  args <- list(
    data = data,
    cutpoints = 0,
    end_of_study = 36,
    prior = c(0.1, 0.1),
    N_mcmc = 10,
    single_arm = FALSE,
    method = "cox",
    alternative = "less"
  )
  res_superiority <- do.call(analyse_data, c(args, h0 = 0))
  res_noninferiority <- do.call(analyse_data, c(args, h0 = log(1.5)))

  expect_true(res_superiority$success < 0.05)
  expect_true(res_noninferiority$success > 0.95)
})
