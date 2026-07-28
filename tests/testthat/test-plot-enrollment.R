test_that("plot_enrollment reproduces constant-rate milestone means", {
  out <- plot_enrollment(
    lambda = 20,
    N_total = 600,
    interim_look = 400,
    end_of_study = 12,
    n_sim = 0,
    time_unit = "months"
  )

  expect_equal(out$milestones$N, c(400, 600))
  expect_equal(out$milestones$projected_time, c(399, 599) / 20)
  expect_equal(out$projection$enrolled[1], 1)
  expect_equal(tail(out$projection$enrolled, 1), 600)
  expect_length(out$simulations, 0)
})

test_that("plot_enrollment supports piecewise enrollment projections", {
  out <- plot_enrollment(
    lambda = c(1, 2),
    lambda_time = 10,
    N_total = 20,
    interim_look = 5,
    n_sim = 0
  )

  expect_equal(out$milestones$projected_time, c(4, 14.5))
  expect_equal(tail(out$projection$enrolled, 1), 20)
})

test_that("plot_enrollment extracts and overrides a retained design", {
  result <- data.frame(N_enrolled = 20)
  attr(result, "enrollment_design") <- new_enrollment_design(
    lambda = 2,
    N_total = 20,
    lambda_time = NULL,
    interim_look = 10,
    end_of_study = 12
  )

  out <- plot_enrollment(result, lambda = 4, n_sim = 0)

  expect_equal(out$design$lambda, 4)
  expect_equal(out$milestones$projected_time, c(9, 19) / 4)
  expect_equal(out$design$end_of_study, 12)
})

test_that("trial results retain evaluated enrollment designs", {
  set.seed(8821)
  result <- survival_adapt(
    hazard_treatment = -log(0.8) / 12,
    hazard_control = NULL,
    N_total = 20,
    lambda = c(2, 4),
    lambda_time = 3,
    interim_look = 10,
    end_of_study = 12,
    alternative = "less",
    h0 = 0.3,
    N_impute = 1,
    N_mcmc = 1,
    method = "bayes-surv",
    return_trace = TRUE
  )

  design <- attr(result, "enrollment_design", exact = TRUE)
  expect_equal(design$lambda, c(2, 4))
  expect_equal(design$lambda_time, 3)
  expect_equal(design$interim_look, 10)
  expect_equal(design$N_total, 20)
  expect_equal(design$end_of_study, 12)

  expect_no_error(plot_enrollment(result, n_sim = 0))
})

test_that("seeded enrollment plots preserve the caller RNG state", {
  set.seed(9013)
  before_kind <- RNGkind()
  before_seed <- .Random.seed

  first <- plot_enrollment(
    lambda = 3,
    N_total = 20,
    n_sim = 2,
    seed = 2026
  )
  expect_identical(RNGkind(), before_kind)
  expect_identical(.Random.seed, before_seed)

  second <- plot_enrollment(
    lambda = 3,
    N_total = 20,
    n_sim = 2,
    seed = 2026
  )
  expect_identical(first$simulations, second$simulations)
})

test_that("plot_enrollment validates its inputs", {
  expect_error(
    plot_enrollment(n_sim = 0),
    "both 'lambda' and 'N_total'"
  )
  expect_error(
    plot_enrollment(lambda = 2, N_total = 20, interim_look = 20, n_sim = 0),
    "less than 'N_total'"
  )
  expect_error(
    plot_enrollment(lambda = 2, N_total = 20, n_sim = -1),
    "non-negative integer"
  )
  expect_error(
    plot_enrollment(lambda = 2, N_total = 20, seed = -1),
    "non-negative integer"
  )
  expect_error(
    plot_enrollment(data.frame(x = 1), n_sim = 0),
    "does not contain a stored"
  )
})
