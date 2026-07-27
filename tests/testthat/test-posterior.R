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
    prior_surv = c(0.1, 0.1),
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
    prior_surv = c(0.1, 0.1),
    N_mcmc = 30,
    single_arm = FALSE
  )
  expect_equal(dim(res), c(30, 3, 2))
})

test_that("vector survival priors broadcast exactly across intervals", {
  data <- data.frame(
    time = c(3, 8, 15, 4, 10, 18),
    event = c(1, 0, 1, 0, 1, 1),
    treatment = c(0, 0, 0, 1, 1, 1)
  )
  prior_matrix <- matrix(rep(c(0.5, 0.25), 3), nrow = 2)

  set.seed(3817)
  from_vector <- posterior(
    data = data,
    cutpoints = c(5, 12),
    prior_surv = c(0.5, 0.25),
    N_mcmc = 50,
    single_arm = FALSE
  )
  set.seed(3817)
  from_matrix <- posterior(
    data = data,
    cutpoints = c(5, 12),
    prior_surv = prior_matrix,
    N_mcmc = 50,
    single_arm = FALSE
  )

  expect_identical(from_matrix, from_vector)
})

test_that("posterior applies the assigned prior to each interval", {
  data <- data.frame(
    time = c(3, 8, 15, 4, 10, 18),
    event = c(1, 0, 1, 0, 1, 1),
    treatment = c(0, 0, 0, 1, 1, 1)
  )
  data_summ <- posterior_sufficient_stats(
    data,
    cutpoints = c(5, 12),
    single_arm = FALSE
  )
  prior_surv <- rbind(
    shape = c(0.5, 2, 5),
    rate = c(0.25, 1, 4)
  )

  set.seed(9017)
  actual <- posterior_from_sufficient_stats(
    data_summ = data_summ,
    prior_surv = prior_surv,
    N_mcmc = 20,
    single_arm = FALSE
  )

  expected <- array(NA_real_, dim = c(20, 3, 2))
  set.seed(9017)
  for (arm_slice in seq_len(2)) {
    treatment_value <- if (arm_slice == 1L) 1 else 0
    arm_summ <- data_summ[data_summ$treatment == treatment_value, ]
    for (j in seq_len(3)) {
      expected[, j, arm_slice] <- rgamma(
        20,
        shape = prior_surv[1, j] + arm_summ$tot_events[j],
        rate = prior_surv[2, j] + arm_summ$tot_time[j]
      )
    }
  }

  expect_identical(actual, expected)
})

test_that("single-interval survival prior matrices are supported", {
  data <- data.frame(
    time = c(2, 4, 3, 5),
    event = c(1, 0, 0, 1),
    treatment = c(0, 0, 1, 1)
  )

  set.seed(749)
  from_vector <- posterior(
    data,
    cutpoints = NULL,
    prior_surv = c(2, 0.5),
    N_mcmc = 20,
    single_arm = FALSE
  )
  set.seed(749)
  from_matrix <- posterior(
    data,
    cutpoints = NULL,
    prior_surv = matrix(c(2, 0.5), nrow = 2),
    N_mcmc = 20,
    single_arm = FALSE
  )

  expect_identical(from_matrix, from_vector)
})

test_that("posterior rejects malformed interval-specific priors", {
  data <- data.frame(
    time = c(3, 8, 15, 4, 10, 18),
    event = c(1, 0, 1, 0, 1, 1),
    treatment = c(0, 0, 0, 1, 1, 1)
  )

  expect_error(
    posterior(
      data,
      cutpoints = c(5, 12),
      prior_surv = matrix(1, nrow = 2, ncol = 2),
      N_mcmc = 10,
      single_arm = FALSE
    ),
    "2 x 3"
  )
  expect_error(
    posterior(
      data,
      cutpoints = c(5, 12),
      prior_surv = rbind(c(1, 1, 1), c(1, 0, 1)),
      N_mcmc = 10,
      single_arm = FALSE
    ),
    "positive finite"
  )
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
    prior_surv = c(0.1, 0.1),
    N_mcmc = 40,
    single_arm = TRUE
  )
  expect_equal(dim(res), c(40, 1, 2))
  expect_true(all(res[,, 1] > 0))
  # Control slice should be all NA for single-arm
  expect_true(all(is.na(res[,, 2])))
})

test_that("statistics-based posterior exactly matches patient-level posterior", {
  # Exercise every shape required by issue #38. Each piecewise example contains
  # exposure in every arm/interval so that this test compares the posterior
  # paths themselves without invoking an empty-interval policy.
  cases <- list(
    "two-arm, single interval" = list(
      data = data.frame(
        time = c(3, 7, 14, 4, 9, 16),
        event = c(1, 0, 1, 0, 1, 0),
        treatment = c(0, 0, 0, 1, 1, 1)
      ),
      cutpoints = NULL,
      single_arm = FALSE
    ),
    "two-arm, piecewise" = list(
      data = data.frame(
        time = c(3, 8, 15, 4, 10, 18),
        event = c(1, 0, 1, 0, 1, 1),
        treatment = c(0, 0, 0, 1, 1, 1)
      ),
      cutpoints = c(5, 12),
      single_arm = FALSE
    ),
    "single-arm, single interval" = list(
      data = data.frame(
        time = c(2, 6, 11, 15),
        event = c(1, 0, 1, 0),
        treatment = rep(1, 4)
      ),
      cutpoints = NULL,
      single_arm = TRUE
    ),
    "single-arm, piecewise" = list(
      data = data.frame(
        time = c(3, 8, 14, 19),
        event = c(1, 0, 1, 0),
        treatment = rep(1, 4)
      ),
      cutpoints = c(5, 12),
      single_arm = TRUE
    )
  )

  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    data_summ <- posterior_sufficient_stats(
      data = case$data,
      cutpoints = case$cutpoints,
      single_arm = case$single_arm
    )

    # No random draws occur while sufficient statistics are calculated.
    # Resetting to the same seed therefore requires exact, not approximate,
    # equality if both posterior entry points preserve draw ordering.
    set.seed(3817)
    from_data <- posterior(
      data = case$data,
      cutpoints = case$cutpoints,
      prior_surv = c(0.5, 0.25),
      N_mcmc = 50,
      single_arm = case$single_arm
    )
    set.seed(3817)
    from_stats <- posterior_from_sufficient_stats(
      data_summ = data_summ,
      prior_surv = c(0.5, 0.25),
      N_mcmc = 50,
      single_arm = case$single_arm
    )

    expect_identical(from_stats, from_data, info = case_name)
  }
})

test_that("posterior sufficient statistics can select analysis rows directly", {
  data <- data.frame(
    time = c(3, 7, 14, 4, 9, 16, 20, 22),
    event = c(1, 0, 1, 0, 1, 0, 1, 0),
    treatment = rep(c(0, 1), each = 4),
    subject_enrolled = c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE)
  )

  from_rows <- posterior_sufficient_stats(
    data = data,
    cutpoints = c(5, 12),
    single_arm = FALSE,
    rows = data$subject_enrolled
  )
  from_subset <- posterior_sufficient_stats(
    data = data[data$subject_enrolled, , drop = FALSE],
    cutpoints = c(5, 12),
    single_arm = FALSE
  )

  expect_identical(from_rows, from_subset)
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
      prior_surv = c(0.1, 0.1),
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
      prior_surv = c(0.1, 0.1),
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
    prior_surv = rbind(
      shape = c(0.5, 1, 5),
      rate = c(0.25, 0.5, 4)
    ),
    N_mcmc = 10000,
    single_arm = FALSE,
    empty_interval = "prior"
  )

  # Treatment interval 3 has no exposure, so with empty_interval = "prior"
  # the posterior mean should be close to the Gamma prior mean.
  expect_equal(mean(res[, 3, 1]), 5 / 4, tolerance = 0.05)
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
      prior_surv = c(0.1, 0.1),
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
      prior_surv = c(0.1, 0.1),
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
      prior_surv = c(0.1, 0.1),
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
      prior_surv = c(0.1, 0.1),
      N_mcmc = 10,
      single_arm = FALSE
    ),
    "No subjects in the control arm"
  )

  # The same all-treatment data is valid for a single-arm analysis.
  res <- posterior(
    data_no_ctrl,
    cutpoints = NULL,
    prior_surv = c(0.1, 0.1),
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
    prior_surv = c(0.1, 0.1),
    N_mcmc = 100,
    single_arm = FALSE
  )
  expect_true(all(res > 0))
})
