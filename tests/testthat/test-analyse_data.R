test_that("analyse_data works with method = 'logrank'", {
  set.seed(6184)
  data <- data.frame(
    time      = rexp(100, 0.02),
    event     = sample(0:1, 100, replace = TRUE),
    treatment = rep(0:1, each = 50)
  )
  res <- analyse_data(
    data         = data,
    cutpoints    = 0,
    end_of_study = 36,
    prior        = c(0.1, 0.1),
    N_mcmc       = 10,
    single_arm   = FALSE,
    method       = "logrank",
    alternative  = "greater",
    h0           = 0
  )
  expect_type(res, "list")
  expect_named(res, c("success", "effect"))
  expect_true(res$success >= 0 && res$success <= 1)
  expect_true(is.na(res$effect))
})

test_that("analyse_data works with method = 'cox'", {
  set.seed(3729)
  data <- data.frame(
    time      = rexp(100, 0.02),
    event     = sample(0:1, 100, replace = TRUE),
    treatment = rep(0:1, each = 50)
  )
  res <- analyse_data(
    data         = data,
    cutpoints    = 0,
    end_of_study = 36,
    prior        = c(0.1, 0.1),
    N_mcmc       = 10,
    single_arm   = FALSE,
    method       = "cox",
    alternative  = "greater",
    h0           = 0
  )
  expect_type(res, "list")
  expect_named(res, c("success", "effect"))
  expect_true(res$success >= 0 && res$success <= 1)
  expect_type(res$effect, "double")
})

test_that("analyse_data works with method = 'bayes' (two-arm)", {
  set.seed(8415)
  data <- data.frame(
    time      = rexp(100, 0.02),
    event     = sample(0:1, 100, replace = TRUE),
    treatment = rep(0:1, each = 50)
  )
  res <- analyse_data(
    data         = data,
    cutpoints    = 0,
    end_of_study = 36,
    prior        = c(0.1, 0.1),
    N_mcmc       = 100,
    single_arm   = FALSE,
    method       = "bayes",
    alternative  = "greater",
    h0           = 0
  )
  expect_type(res, "list")
  expect_named(res, c("success", "effect"))
  expect_true(res$success >= 0 && res$success <= 1)
  expect_length(res$effect, 100)
})

test_that("analyse_data works with method = 'bayes' and alternative = 'less'", {
  set.seed(5091)
  data <- data.frame(
    time      = rexp(100, 0.02),
    event     = sample(0:1, 100, replace = TRUE),
    treatment = rep(0:1, each = 50)
  )
  res_greater <- analyse_data(
    data         = data,
    cutpoints    = 0,
    end_of_study = 36,
    prior        = c(0.1, 0.1),
    N_mcmc       = 500,
    single_arm   = FALSE,
    method       = "bayes",
    alternative  = "greater",
    h0           = 0
  )
  res_less <- analyse_data(
    data         = data,
    cutpoints    = 0,
    end_of_study = 36,
    prior        = c(0.1, 0.1),
    N_mcmc       = 500,
    single_arm   = FALSE,
    method       = "bayes",
    alternative  = "less",
    h0           = 0
  )
  # P(effect > 0) + P(effect < 0) should be approximately 1
  # (exact only if no draws == 0)
  expect_equal(res_greater$success + res_less$success, 1, tolerance = 0.05)
})

test_that("analyse_data works with method = 'chisq'", {
  set.seed(2647)
  data <- data.frame(
    time      = rexp(100, 0.02),
    event     = sample(0:1, 100, replace = TRUE, prob = c(0.3, 0.7)),
    treatment = rep(0:1, each = 50)
  )
  res <- analyse_data(
    data         = data,
    cutpoints    = 0,
    end_of_study = 36,
    prior        = c(0.1, 0.1),
    N_mcmc       = 10,
    single_arm   = FALSE,
    method       = "chisq",
    alternative  = "greater",
    h0           = 0
  )
  expect_type(res, "list")
  expect_named(res, c("success", "effect"))
  expect_true(res$success >= 0 && res$success <= 1)
  expect_type(res$effect, "double")
})

test_that("analyse_data works with method = 'bayes' (single-arm)", {
  set.seed(9473)
  data <- data.frame(
    time      = rexp(50, 0.02),
    event     = sample(0:1, 50, replace = TRUE),
    treatment = rep(1, 50)
  )
  res <- analyse_data(
    data         = data,
    cutpoints    = 0,
    end_of_study = 36,
    prior        = c(0.1, 0.1),
    N_mcmc       = 100,
    single_arm   = TRUE,
    method       = "bayes",
    alternative  = "greater",
    h0           = 0
  )
  expect_type(res, "list")
  expect_true(res$success >= 0 && res$success <= 1)
  expect_length(res$effect, 100)
})

test_that("analyse_data works with piecewise cutpoints", {
  set.seed(4318)
  data <- data.frame(
    time      = rexp(100, 0.02),
    event     = sample(0:1, 100, replace = TRUE),
    treatment = rep(0:1, each = 50)
  )
  res <- analyse_data(
    data         = data,
    cutpoints    = c(0, 12),
    end_of_study = 36,
    prior        = c(0.1, 0.1),
    N_mcmc       = 100,
    single_arm   = FALSE,
    method       = "bayes",
    alternative  = "greater",
    h0           = 0
  )
  expect_type(res, "list")
  expect_true(res$success >= 0 && res$success <= 1)
})
