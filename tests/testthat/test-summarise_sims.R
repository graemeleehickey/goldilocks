test_that("summarise_sims works with a single data frame", {
  df <- data.frame(
    stop_futility = c(FALSE, FALSE, TRUE, FALSE),
    post_prob_ha = c(0.98, 0.80, 0.30, 0.99),
    prob_threshold = c(0.95, 0.95, 0.95, 0.95),
    stop_expected_success = c(TRUE, FALSE, FALSE, TRUE),
    N_enrolled = c(200, 400, 150, 300)
  )
  out <- summarise_sims(df)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 1)
  expect_true("power" %in% names(out))
  expect_true("mean_N" %in% names(out))
  expect_true("sd_N" %in% names(out))
  expect_identical(out$n_analyzed, 4L)
  expect_identical(out$n_used, 4L)
  expect_true(is.na(out$n_requested))
  expect_true(is.na(out$n_failed))
  expect_true(is.na(out$failure_rate))
  expect_equal(attr(out, "mc_conf_level", exact = TRUE), 0.95)
})

test_that("summarise_sims works with a named list of data frames", {
  make_df <- function(n) {
    data.frame(
      stop_futility = sample(c(TRUE, FALSE), n, replace = TRUE),
      post_prob_ha = runif(n),
      prob_threshold = rep(0.95, n),
      stop_expected_success = sample(c(TRUE, FALSE), n, replace = TRUE),
      N_enrolled = sample(100:500, n, replace = TRUE)
    )
  }
  set.seed(3059)
  data_list <- list(null = make_df(10), alt = make_df(10))
  out <- summarise_sims(data_list)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 2)
  expect_setequal(out$scenario, c("null", "alt"))
})

test_that("summarise_sims labels unnamed list scenarios by position", {
  make_df <- function(n) {
    data.frame(
      stop_futility = rep(FALSE, n),
      post_prob_ha = rep(0.99, n),
      prob_threshold = rep(0.95, n),
      stop_expected_success = rep(FALSE, n),
      N_enrolled = seq_len(n)
    )
  }
  out <- summarise_sims(list(make_df(3), make_df(4)))
  expect_s3_class(out, "data.frame")
  expect_equal(out$scenario, c("1", "2"))
})

test_that("summarise_sims computes power correctly", {
  # All trials succeed: no futility, posterior > threshold
  df <- data.frame(
    stop_futility = rep(FALSE, 5),
    post_prob_ha = rep(0.99, 5),
    prob_threshold = rep(0.95, 5),
    stop_expected_success = rep(TRUE, 5),
    N_enrolled = rep(200, 5)
  )
  out <- summarise_sims(df)
  expect_equal(out$power, 1)
  expect_lt(out$power_mc_lower, 1)
  expect_equal(out$power_mc_upper, 1)
})

test_that("summarise_sims computes zero power when all futile", {
  df <- data.frame(
    stop_futility = rep(TRUE, 5),
    post_prob_ha = rep(0.10, 5),
    prob_threshold = rep(0.95, 5),
    stop_expected_success = rep(FALSE, 5),
    N_enrolled = rep(100, 5)
  )
  out <- summarise_sims(df)
  expect_equal(out$power, 0)
  expect_equal(out$power_mc_lower, 0)
  expect_gt(out$power_mc_upper, 0)
})

test_that("Wilson intervals handle zero and one observed outcome", {
  uncertainty <- binomial_mc_uncertainty(
    estimate = c(0, 0.1, 1),
    denominator = c(10, 10, 10)
  )
  z <- stats::qnorm(0.975)
  boundary_width <- z^2 / (10 + z^2)

  expect_equal(uncertainty$lower[1], 0)
  expect_equal(uncertainty$upper[1], boundary_width)
  expect_gt(uncertainty$lower[2], 0)
  expect_lt(uncertainty$upper[2], 1)
  expect_equal(uncertainty$lower[3], 1 - boundary_width)
  expect_equal(uncertainty$upper[3], 1)
})

test_that("mean sample size uses a t Monte Carlo interval", {
  df <- data.frame(
    stop_futility = rep(FALSE, 4),
    post_prob_ha = rep(0.99, 4),
    prob_threshold = rep(0.95, 4),
    stop_expected_success = rep(TRUE, 4),
    N_enrolled = c(100, 150, 200, 250)
  )
  out <- summarise_sims(df)
  expected_mcse <- stats::sd(df$N_enrolled) / sqrt(4)
  half_width <- stats::qt(0.975, df = 3) * expected_mcse

  expect_equal(out$mean_N_mcse, expected_mcse)
  expect_equal(out$mean_N_mc_lower, mean(df$N_enrolled) - half_width)
  expect_equal(out$mean_N_mc_upper, mean(df$N_enrolled) + half_width)
})

test_that("maximum sample size is calculated from trial-level indicators", {
  df <- data.frame(
    stop_futility = c(FALSE, FALSE, TRUE, FALSE),
    post_prob_ha = c(0.99, 0.80, 0.30, 0.99),
    prob_threshold = rep(0.95, 4),
    stop_expected_success = c(TRUE, FALSE, FALSE, TRUE),
    N_enrolled = c(200, 400, 150, 300)
  )
  out <- summarise_sims(df)

  expect_equal(out$stop_success, 0.5)
  expect_equal(out$stop_futility, 0.25)
  expect_equal(out$stop_max_N, 0.25)
  expect_equal(
    out$stop_success + out$stop_futility + out$stop_max_N,
    1
  )
})

test_that("summarise_sims accepts a complete sim_trials result", {
  simulations <- sim_trials(
    hazard_treatment = 0.02,
    hazard_control = 0.03,
    N_total = 20,
    lambda = 10,
    end_of_study = 12,
    N_trials = 3,
    method = "logrank",
    backend = "sequential",
    seed = 7801
  )

  direct <- summarise_sims(simulations)
  legacy <- summarise_sims(simulations$sims)
  metric_columns <- c(
    "power",
    "stop_success",
    "stop_futility",
    "stop_max_N",
    "mean_N",
    "sd_N",
    "stop_and_fail"
  )

  for (column in metric_columns) {
    expect_equal(direct[[column]], legacy[[column]])
  }
  expect_identical(direct$backend, "sequential")
  expect_equal(direct$seed, 7801)
  expect_identical(direct$n_requested, 3L)
  expect_identical(direct$n_analyzed, 3L)
  expect_identical(direct$n_failed, 0L)
  expect_identical(direct$n_used, 3L)
  expect_equal(direct$failure_rate, 0)
  expect_true(all(c(
    "power_mcse",
    "power_mc_lower",
    "power_mc_upper",
    "mean_N_mcse",
    "mean_N_mc_lower",
    "mean_N_mc_upper"
  ) %in% names(direct)))
  expect_identical(
    attr(direct, "enrollment_design", exact = TRUE),
    attr(simulations, "enrollment_design", exact = TRUE)
  )
  expect_identical(
    attr(direct, "decision_design", exact = TRUE),
    attr(simulations, "decision_design", exact = TRUE)
  )
  expect_named(attr(direct, "simulation_metadata", exact = TRUE), "1")
})

test_that("summarise_sims retains named full-result scenarios", {
  make_result <- function(seed) {
    sim_trials(
      hazard_treatment = 0.02,
      hazard_control = 0.03,
      N_total = 20,
      lambda = 10,
      end_of_study = 12,
      N_trials = 2,
      method = "logrank",
      backend = "sequential",
      seed = seed
    )
  }

  out <- summarise_sims(list(
    null = make_result(7802),
    alternative = make_result(7803)
  ))

  expect_setequal(out$scenario, c("null", "alternative"))
  expect_equal(
    out$seed[match(c("null", "alternative"), out$scenario)],
    c(7802, 7803)
  )
  expect_identical(out$n_requested, c(2L, 2L))
  expect_named(
    attr(out, "simulation_metadata", exact = TRUE),
    c("null", "alternative")
  )
})

test_that("summarise_sims retains grouped identifiers", {
  data <- data.frame(
    cohort = rep(c("A", "B"), each = 2),
    scenario = rep(c("null", "target"), 2),
    stop_futility = c(TRUE, FALSE, FALSE, FALSE),
    post_prob_ha = c(0.1, 0.99, 0.2, 0.98),
    prob_threshold = rep(0.95, 4),
    stop_expected_success = c(FALSE, TRUE, FALSE, TRUE),
    N_enrolled = c(10, 20, 15, 20)
  ) |>
    dplyr::group_by(cohort)

  out <- summarise_sims(data)

  expect_true(all(c("cohort", "scenario") %in% names(out)))
  expect_equal(nrow(out), 4L)
  expect_setequal(out$cohort, c("A", "B"))
  expect_setequal(out$scenario, c("null", "target"))
  expect_identical(out$n_used, rep(1L, 4))
  expect_true(all(is.finite(out$power_mc_lower)))
  expect_true(all(is.finite(out$power_mc_upper)))
  expect_true(all(is.na(out$mean_N_mcse)))
})

test_that("summarise_sims reports failures separately from estimand denominators", {
  simulations <- sim_trials(
    hazard_treatment = 0.02,
    hazard_control = 0.03,
    N_total = 20,
    lambda = 10,
    end_of_study = 12,
    N_trials = 3,
    method = "logrank",
    backend = "sequential",
    seed = 7804
  )
  simulations$sims <- simulations$sims[1:2, , drop = FALSE]
  simulations$sims$stop_futility <- c(FALSE, TRUE)
  simulations$sims$stop_expected_success <- c(TRUE, FALSE)
  simulations$sims$post_prob_ha <- c(0.99, 0.10)
  simulations$sims$prob_threshold <- c(0.95, 0.95)
  simulations$failures <- data.frame(
    trial = 3:4,
    message = c("failure one", "failure two")
  )
  parallel_metadata <- attr(simulations, "parallel_metadata", exact = TRUE)
  parallel_metadata$tasks <- 4L
  attr(simulations, "parallel_metadata") <- parallel_metadata

  out <- summarise_sims(simulations)
  successful <- with(
    simulations$sims,
    !stop_futility & post_prob_ha > prob_threshold
  )

  expect_identical(out$n_requested, 4L)
  expect_identical(out$n_analyzed, 2L)
  expect_identical(out$n_failed, 2L)
  expect_identical(out$n_used, 2L)
  expect_equal(out$failure_rate, 0.5)
  expect_equal(out$power, sum(successful) / 2)
  expect_equal(
    out$stop_success,
    sum(simulations$sims$stop_expected_success) / 2
  )
  expect_equal(
    attr(out, "simulation_metadata", exact = TRUE)[[1]]$failures,
    simulations$failures
  )
})

test_that("summarise_sims warns when requested precision is not achieved", {
  df <- data.frame(
    stop_futility = c(FALSE, FALSE, TRUE, TRUE),
    post_prob_ha = c(0.99, 0.99, 0.10, 0.10),
    prob_threshold = rep(0.95, 4),
    stop_expected_success = c(TRUE, TRUE, FALSE, FALSE),
    N_enrolled = c(100, 120, 140, 160)
  )

  expect_warning(
    summarise_sims(df, max_mcse = c(power = 0.1)),
    "scenario=1, power: MCSE 0.2500 exceeds 0.1000"
  )
  expect_no_warning(
    summarise_sims(df, max_mcse = c(power = 0.3, mean_N = 20))
  )
})

test_that("summarise_sims validates precision targets", {
  df <- data.frame(
    stop_futility = FALSE,
    post_prob_ha = 0.99,
    prob_threshold = 0.95,
    stop_expected_success = TRUE,
    N_enrolled = 100
  )

  expect_error(summarise_sims(df, max_mcse = 0.1), "named numeric")
  expect_error(
    summarise_sims(df, max_mcse = c(unknown = 0.1)),
    "Unsupported.*unknown"
  )
  expect_error(
    summarise_sims(df, max_mcse = c(power = 0)),
    "positive"
  )
})

test_that("summarise_sims rejects invalid list-like inputs clearly", {
  expect_error(
    summarise_sims(list(valid = data.frame(), invalid = list(value = 1))),
    "invalid list element.*position.*2"
  )
  expect_error(
    summarise_sims(list()),
    "non-empty list"
  )
  expect_error(
    summarise_sims(list(sims = data.frame(), call = "not a call")),
    "invalid list element"
  )
})
