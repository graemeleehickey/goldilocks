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

test_that("events at cutpoints use open-left closed-right intervals", {
  epsilon <- 1e-6
  data <- data.frame(
    time = c(
      5 - epsilon,
      5,
      5 + epsilon,
      12 - epsilon,
      12,
      12 + epsilon
    ),
    event = 1,
    treatment = 1
  )

  stats <- posterior_sufficient_stats(
    data,
    cutpoints = c(5, 12),
    single_arm = TRUE
  )

  expect_equal(stats$tot_events, c(2, 3, 1))
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
      single_arm = FALSE,
      empty_interval = "propagate"
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
      single_arm = FALSE,
      empty_interval = "propagate"
    )
  )

  trt_int2 <- mean(res[, 2, 1])
  trt_int3 <- mean(res[, 3, 1])
  ctrl_int3 <- mean(res[, 3, 2])

  # Back-filled treatment interval is closer to its own group than to control.
  expect_lt(abs(trt_int3 - trt_int2), abs(trt_int3 - ctrl_int3))
})

test_that("posterior leaves zero-exposure intervals prior-driven by default", {
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
    single_arm = FALSE
  )

  # Treatment interval 3 has no exposure, so with empty_interval = "prior"
  # the posterior mean should be close to the Gamma prior mean.
  expect_equal(mean(res[, 3, 1]), 5 / 4, tolerance = 0.05)
})

test_that("prior-only handling covers empty first, middle, and final intervals", {
  prior_surv <- rbind(
    shape = c(2, 6, 12),
    rate = c(4, 3, 4)
  )
  data_summ <- data.frame(
    treatment = c(0, 0, 0, 1, 1, 1),
    interval = factor(rep(1:3, 2), levels = 1:3),
    n = c(2, 0, 2, 0, 2, 0),
    tot_time = c(5, 0, 8, 0, 6, 0),
    tot_events = c(1, 0, 1, 0, 1, 0)
  )

  set.seed(621)
  post <- posterior_from_sufficient_stats(
    data_summ = data_summ,
    prior_surv = prior_surv,
    N_mcmc = 10000,
    single_arm = FALSE
  )

  # Control interval 2 and treatment intervals 1 and 3 are empty. Their
  # posterior means remain the corresponding interval-specific prior means.
  expect_equal(mean(post[, 2, 2]), 6 / 3, tolerance = 0.05)
  expect_equal(mean(post[, 1, 1]), 2 / 4, tolerance = 0.03)
  expect_equal(mean(post[, 3, 1]), 12 / 4, tolerance = 0.07)
})

test_that("single-arm prior-only handling covers every empty interval position", {
  prior_surv <- rbind(
    shape = c(2, 6, 12),
    rate = c(4, 3, 4)
  )

  for (empty_interval in 1:3) {
    data_summ <- data.frame(
      treatment = rep(1, 3),
      interval = factor(1:3, levels = 1:3),
      n = rep(2, 3),
      tot_time = rep(5, 3),
      tot_events = rep(1, 3)
    )
    data_summ[empty_interval, c("n", "tot_time", "tot_events")] <- 0

    set.seed(630 + empty_interval)
    post <- posterior_from_sufficient_stats(
      data_summ = data_summ,
      prior_surv = prior_surv,
      N_mcmc = 8000,
      single_arm = TRUE
    )

    expect_equal(
      mean(post[, empty_interval, 1]),
      unname(
        prior_surv[1, empty_interval] / prior_surv[2, empty_interval]
      ),
      tolerance = 0.07,
      info = paste("empty interval", empty_interval)
    )
  }
})

test_that("prior-only empty-interval results are invariant to row ordering", {
  data_summ <- data.frame(
    treatment = c(0, 0, 0, 1, 1, 1),
    interval = factor(rep(1:3, 2), levels = 1:3),
    n = c(2, 0, 2, 0, 2, 0),
    tot_time = c(5, 0, 8, 0, 6, 0),
    tot_events = c(1, 0, 1, 0, 1, 0)
  )
  shuffled <- data_summ[c(6, 2, 4, 1, 5, 3), ]

  set.seed(622)
  ordered_post <- posterior_from_sufficient_stats(
    data_summ,
    prior_surv = c(0.5, 0.25),
    N_mcmc = 100,
    single_arm = FALSE
  )
  set.seed(622)
  shuffled_post <- posterior_from_sufficient_stats(
    shuffled,
    prior_surv = c(0.5, 0.25),
    N_mcmc = 100,
    single_arm = FALSE
  )

  expect_identical(shuffled_post, ordered_post)
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

test_that("posterior applies independent arm-specific priors", {
  data <- data.frame(
    time = c(3, 8, 15, 4, 10, 18),
    event = c(1, 0, 1, 0, 1, 1),
    treatment = c(0, 0, 0, 1, 1, 1)
  )
  cutpoints <- c(5, 12)
  data_summ <- posterior_sufficient_stats(
    data,
    cutpoints = cutpoints,
    single_arm = FALSE
  )
  prior_surv <- list(
    treatment = rbind(
      shape = c(5, 6, 7),
      rate = c(8, 9, 10)
    ),
    control = rbind(
      shape = c(0.5, 1, 2),
      rate = c(0.25, 0.5, 1)
    )
  )

  set.seed(8891)
  actual <- posterior_from_sufficient_stats(
    data_summ = data_summ,
    prior_surv = prior_surv,
    N_mcmc = 20,
    single_arm = FALSE
  )

  expected <- array(NA_real_, dim = c(20, 3, 2))
  set.seed(8891)
  for (arm_slice in seq_len(2)) {
    arm <- if (arm_slice == 1L) "treatment" else "control"
    treatment_value <- if (arm == "treatment") 1 else 0
    arm_summ <- data_summ[data_summ$treatment == treatment_value, ]
    for (j in seq_len(3)) {
      expected[, j, arm_slice] <- rgamma(
        20,
        shape = prior_surv[[arm]]["shape", j] + arm_summ$tot_events[j],
        rate = prior_surv[[arm]]["rate", j] + arm_summ$tot_time[j]
      )
    }
  }

  expect_identical(actual, expected)
})

test_that("equal arm-specific priors preserve legacy seeded draws", {
  data <- data.frame(
    time = c(3, 8, 15, 4, 10, 18),
    event = c(1, 0, 1, 0, 1, 1),
    treatment = c(0, 0, 0, 1, 1, 1)
  )
  shared <- rbind(
    shape = c(0.5, 2, 5),
    rate = c(0.25, 1, 4)
  )

  set.seed(8892)
  legacy <- posterior(
    data,
    cutpoints = c(5, 12),
    prior_surv = shared,
    N_mcmc = 50,
    single_arm = FALSE
  )
  set.seed(8892)
  arm_specific <- posterior(
    data,
    cutpoints = c(5, 12),
    prior_surv = list(treatment = shared, control = shared),
    N_mcmc = 50,
    single_arm = FALSE
  )

  expect_identical(arm_specific, legacy)
})

test_that("changing one arm's prior does not alter the other arm's draws", {
  data <- data.frame(
    time = c(3, 8, 15, 4, 10, 18),
    event = c(1, 0, 1, 0, 1, 1),
    treatment = c(0, 0, 0, 1, 1, 1)
  )
  shared_treatment <- c(2, 4)

  set.seed(8894)
  weak_control <- posterior(
    data,
    cutpoints = c(5, 12),
    prior_surv = list(
      control = c(0.1, 0.1),
      treatment = shared_treatment
    ),
    N_mcmc = 50,
    single_arm = FALSE
  )
  set.seed(8894)
  strong_control <- posterior(
    data,
    cutpoints = c(5, 12),
    prior_surv = list(
      control = c(100, 100),
      treatment = shared_treatment
    ),
    N_mcmc = 50,
    single_arm = FALSE
  )

  expect_identical(strong_control[,, 1], weak_control[,, 1])
  expect_false(identical(strong_control[,, 2], weak_control[,, 2]))
})

test_that("arm-specific prior names determine arms rather than list order", {
  control <- rbind(shape = c(1, 2), rate = c(3, 4))
  treatment <- rbind(shape = c(5, 6), rate = c(7, 8))

  normalized <- goldilocks:::normalize_gamma_prior(
    list(treatment = treatment, control = control),
    n_intervals = 2,
    single_arm = FALSE
  )

  expect_identical(dim(normalized), c(2L, 2L, 2L))
  expect_identical(dimnames(normalized)$arm, c("control", "treatment"))
  expect_equal(unname(normalized[,, "control"]), unname(control))
  expect_equal(unname(normalized[,, "treatment"]), unname(treatment))
})

test_that("arm-specific survival priors require exact complete arm names", {
  expect_error(
    goldilocks:::normalize_gamma_prior(
      list(c(1, 1), c(2, 2)),
      n_intervals = 1,
      single_arm = FALSE
    ),
    "exactly: control, treatment"
  )
  expect_error(
    goldilocks:::normalize_gamma_prior(
      list(control = c(1, 1)),
      n_intervals = 1,
      single_arm = FALSE
    ),
    "exactly: control, treatment"
  )
  expect_error(
    goldilocks:::normalize_gamma_prior(
      list(control = c(1, 1), treatment = c(2, 0)),
      n_intervals = 1,
      single_arm = FALSE
    ),
    "prior_surv\\$treatment.*positive finite"
  )
  expect_error(
    goldilocks:::normalize_gamma_prior(
      list(control = c(1, 1), treatment = c(2, 2)),
      n_intervals = 1,
      single_arm = TRUE
    ),
    "exactly: treatment"
  )

  single_arm <- goldilocks:::normalize_gamma_prior(
    list(treatment = c(2, 4)),
    n_intervals = 2,
    single_arm = TRUE
  )
  expect_identical(dimnames(single_arm)$arm, "treatment")
  expect_equal(unname(single_arm["shape", , "treatment"]), c(2, 2))
  expect_equal(unname(single_arm["rate", , "treatment"]), c(4, 4))
})

test_that("Gamma diagnostics expose resolved conjugate updates", {
  data <- data.frame(
    time = c(3, 8, 15, 4, 10, 18),
    event = c(1, 0, 1, 0, 1, 1),
    treatment = c(0, 0, 0, 1, 1, 1)
  )
  cutpoints <- c(5, 12)
  data_summ <- posterior_sufficient_stats(
    data,
    cutpoints = cutpoints,
    single_arm = FALSE
  )
  prior_surv <- list(
    control = c(2, 4),
    treatment = c(5, 10)
  )

  diagnostics <- goldilocks:::gamma_posterior_diagnostics(
    data_summ = data_summ,
    prior_surv = prior_surv,
    cutpoints = cutpoints,
    end_of_study = 20,
    single_arm = FALSE,
    empty_interval = "prior"
  )

  expect_identical(diagnostics$arm, rep(c("control", "treatment"), each = 3))
  expect_equal(diagnostics$interval_start, rep(c(0, 5, 12), 2))
  expect_equal(diagnostics$interval_end, rep(c(5, 12, 20), 2))
  expect_equal(
    diagnostics$posterior_shape,
    diagnostics$prior_shape + diagnostics$observed_events
  )
  expect_equal(
    diagnostics$posterior_rate,
    diagnostics$prior_rate + diagnostics$observed_exposure
  )
  expect_equal(
    diagnostics$posterior_mean_hazard,
    diagnostics$posterior_shape / diagnostics$posterior_rate
  )
  expect_false(any(c("draw", "draws") %in% names(diagnostics)))
})

test_that("posterior diagnostics distinguish observed and propagated data", {
  data <- data.frame(
    time = c(2, 25, 3, 4),
    event = c(1, 1, 1, 0),
    treatment = c(0, 0, 1, 1)
  )
  data_summ <- posterior_sufficient_stats(
    data,
    cutpoints = c(10, 20),
    single_arm = FALSE
  )

  diagnostics <- expect_no_warning(
    goldilocks:::gamma_posterior_diagnostics(
      data_summ = data_summ,
      prior_surv = c(0.1, 0.1),
      cutpoints = c(10, 20),
      end_of_study = 30,
      single_arm = FALSE,
      empty_interval = "propagate"
    )
  )
  treatment <- diagnostics[diagnostics$arm == "treatment", ]

  expect_identical(treatment$empty_interval, c(FALSE, TRUE, TRUE))
  expect_identical(
    treatment$empty_interval_policy,
    c("observed", "propagate", "propagate")
  )
  expect_equal(treatment$observed_exposure[2:3], c(0, 0))
  expect_equal(
    treatment$effective_exposure[2:3],
    rep(treatment$effective_exposure[1], 2)
  )
  expect_equal(
    treatment$effective_events[2:3],
    rep(treatment$effective_events[1], 2)
  )
})
