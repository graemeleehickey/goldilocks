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
  data <- data.frame(
    time = rexp(100, 0.02),
    event = sample(0:1, 100, replace = TRUE, prob = c(0.3, 0.7)),
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
