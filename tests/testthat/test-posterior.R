test_that("posterior returns correct dimensions (two-arm, single interval)", {
  set.seed(3821)
  data <- data.frame(
    time = rexp(80, 0.02),
    event = sample(0:1, 80, replace = TRUE),
    treatment = rep(0:1, each = 40)
  )
  res <- posterior(
    data,
    cutpoints = NULL,
    prior = c(0.1, 0.1),
    N_mcmc = 50,
    single_arm = FALSE
  )
  expect_true(is.array(res))
  expect_equal(dim(res), c(50, 1, 2))
  expect_true(all(res[,, 1] > 0))
  expect_true(all(res[,, 2] > 0))
})

test_that("posterior returns correct dimensions (two-arm, piecewise)", {
  set.seed(6194)
  data <- data.frame(
    time = rexp(100, 0.02),
    event = sample(0:1, 100, replace = TRUE),
    treatment = rep(0:1, each = 50)
  )
  res <- posterior(
    data,
    cutpoints = c(12, 24),
    prior = c(0.1, 0.1),
    N_mcmc = 30,
    single_arm = FALSE
  )
  expect_equal(dim(res), c(30, 3, 2))
})

test_that("posterior returns correct dimensions (single-arm)", {
  set.seed(7452)
  data <- data.frame(
    time = rexp(50, 0.02),
    event = sample(0:1, 50, replace = TRUE),
    treatment = rep(1, 50)
  )
  res <- posterior(
    data,
    cutpoints = NULL,
    prior = c(0.1, 0.1),
    N_mcmc = 40,
    single_arm = TRUE
  )
  expect_equal(dim(res), c(40, 1, 2))
  expect_true(all(res[,, 1] > 0))
  # Control slice should be all NA for single-arm
  expect_true(all(is.na(res[,, 2])))
})

test_that("posterior warns on zero-exposure interval", {
  set.seed(2913)
  # Treatment event times are all < 10, so treatment interval 2 (>= 10) has
  # zero exposure. Control contributes exposure in interval 2, so only the
  # expected treatment warning is emitted.
  data <- data.frame(
    time = c(runif(30, 12, 18), runif(30, 0, 9)),
    event = sample(0:1, 60, replace = TRUE),
    treatment = rep(0:1, each = 30)
  )
  expect_warning(
    posterior(
      data,
      cutpoints = 10,
      prior = c(0.1, 0.1),
      N_mcmc = 10,
      single_arm = FALSE
    ),
    "zero subjects"
  )
})

test_that("posterior propagates zero-exposure intervals within treatment group", {
  # Only the control group has follow-up reaching interval 3; the treatment
  # group's interval 3 is empty. The empty treatment interval must be
  # back-filled from the treatment group's own data, not from control data.
  data <- data.frame(
    time = c(2, 25, 3, 4),
    event = c(1, 1, 1, 0),
    treatment = c(0, 0, 1, 1)
  )
  set.seed(771)
  res <- suppressWarnings(
    posterior(
      data,
      cutpoints = c(10, 20),
      prior = c(0.1, 0.1),
      N_mcmc = 4000,
      single_arm = FALSE
    )
  )

  trt_int2 <- mean(res[, 2, 1])
  trt_int3 <- mean(res[, 3, 1])
  ctrl_int3 <- mean(res[, 3, 2])

  # Back-filled treatment interval is closer to its own group than to control.
  expect_lt(abs(trt_int3 - trt_int2), abs(trt_int3 - ctrl_int3))
})

test_that("posterior can leave zero-exposure intervals prior-driven", {
  data <- data.frame(
    time = c(2, 25, 3, 4),
    event = c(1, 1, 1, 0),
    treatment = c(0, 0, 1, 1)
  )

  set.seed(1429)
  res <- posterior(
    data,
    cutpoints = c(10, 20),
    prior = c(2, 0.5),
    N_mcmc = 10000,
    single_arm = FALSE,
    empty_interval = "prior"
  )

  # Treatment interval 3 has no exposure, so with empty_interval = "prior"
  # the posterior mean should be close to the Gamma prior mean.
  expect_equal(mean(res[, 3, 1]), 2 / 0.5, tolerance = 0.15)
})

test_that("posterior can error on zero-exposure intervals", {
  data <- data.frame(
    time = c(2, 25, 3, 4),
    event = c(1, 1, 1, 0),
    treatment = c(0, 0, 1, 1)
  )

  expect_error(
    posterior(
      data,
      cutpoints = c(10, 20),
      prior = c(0.1, 0.1),
      N_mcmc = 10,
      single_arm = FALSE,
      empty_interval = "error"
    ),
    "zero subjects"
  )
})

test_that("posterior errors when the treatment group has no subjects (factor)", {
  # Two-arm analysis but the treatment group has no enrolled subjects in this
  # look, encoded as a factor with both levels present.
  data <- data.frame(
    time = c(2, 3, 4, 5),
    event = c(1, 1, 0, 1),
    treatment = factor(c(0, 0, 0, 0), levels = c(0, 1))
  )
  expect_error(
    posterior(
      data,
      cutpoints = 10,
      prior = c(0.1, 0.1),
      N_mcmc = 10,
      single_arm = FALSE
    ),
    "No subjects in the treatment arm"
  )
})

test_that("posterior errors on single-arm interim data with numeric treatment", {
  # This mirrors what survival_adapt() actually passes: `treatment` is a numeric
  # column, and at a small interim look an entire treatment group can be absent.
  # Previously the absent group produced no summary row and posterior() silently
  # returned an all-NA slice via rgamma(N_mcmc, ..., NA). It must now error.

  # Treatment group absent (all control), two-arm analysis.
  data_no_trt <- data.frame(
    time = c(0.5, 1.2),
    event = c(0, 0),
    treatment = c(0, 0)
  )
  expect_error(
    posterior(
      data_no_trt,
      cutpoints = NULL,
      prior = c(0.1, 0.1),
      N_mcmc = 10,
      single_arm = FALSE
    ),
    "No subjects in the treatment arm"
  )

  # Control group absent (all treatment), two-arm analysis.
  data_no_ctrl <- data.frame(
    time = c(0.5, 1.2),
    event = c(0, 1),
    treatment = c(1, 1)
  )
  expect_error(
    posterior(
      data_no_ctrl,
      cutpoints = NULL,
      prior = c(0.1, 0.1),
      N_mcmc = 10,
      single_arm = FALSE
    ),
    "No subjects in the control arm"
  )

  # The same all-treatment data is valid for a single-arm analysis.
  res <- posterior(
    data_no_ctrl,
    cutpoints = NULL,
    prior = c(0.1, 0.1),
    N_mcmc = 10,
    single_arm = TRUE
  )
  expect_true(all(res[,, 1] > 0))
  expect_true(all(is.na(res[,, 2])))
})

test_that("posterior returns positive draws", {
  set.seed(5087)
  data <- data.frame(
    time = rexp(100, 0.03),
    event = sample(0:1, 100, replace = TRUE),
    treatment = rep(0:1, each = 50)
  )
  res <- posterior(
    data,
    cutpoints = 12,
    prior = c(0.1, 0.1),
    N_mcmc = 100,
    single_arm = FALSE
  )
  expect_true(all(res > 0))
})
