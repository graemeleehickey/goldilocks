test_that("enrollment returns the requested continuous-time schedule", {
  set.seed(2198)
  out <- enrollment(lambda = 5, N_total = 100)

  expect_length(out, 100)
  expect_equal(out[1], 0)
  expect_false(is.unsorted(out))
  expect_true(all(diff(out) > 0))
  expect_true(any(out[-1] != floor(out[-1])))
})

test_that("enrollment handles a one-patient trial", {
  expect_equal(enrollment(lambda = 5, N_total = 1), 0)
})

test_that("constant-rate gaps have the exponential-process distribution", {
  set.seed(7413)
  rate <- 4
  gaps <- diff(enrollment(lambda = rate, N_total = 50000))

  expect_lt(abs(mean(gaps) - 1 / rate), 0.004)
  expect_lt(abs(sd(gaps) - 1 / rate), 0.004)
  expect_lt(abs(mean(gaps <= log(2) / rate) - 0.5), 0.01)
})

test_that("fractional internal knots follow the integrated intensity", {
  # For the first arrival after the patient fixed at zero,
  # P(T <= t) = 1 - exp(-Lambda(t)). The knot at 0.75 deliberately does not
  # align with a unit-time boundary.
  set.seed(4821)
  arrivals <- replicate(
    20000,
    enrollment(
      lambda = c(0.5, 2),
      N_total = 2,
      lambda_time = 0.75
    )[2]
  )

  expected_before_knot <- 1 - exp(-0.5 * 0.5)
  expected_after_knot <- 1 - exp(-(0.5 * 0.75 + 2 * (1 - 0.75)))

  expect_lt(abs(mean(arrivals <= 0.5) - expected_before_knot), 0.01)
  expect_lt(abs(mean(arrivals <= 1) - expected_after_knot), 0.01)
})

test_that("enrollment supports multiple internal knots", {
  set.seed(6739)
  out <- enrollment(
    lambda = c(0.3, 0.7, 1.2),
    N_total = 80,
    lambda_time = c(10, 20)
  )

  expect_length(out, 80)
  expect_equal(out[1], 0)
  expect_false(is.unsorted(out))
})

test_that("enrollment validates rates and sample size", {
  expect_error(
    enrollment(lambda = -1, N_total = 50),
    "finite positive"
  )
  expect_error(
    enrollment(lambda = 0, N_total = 50),
    "finite positive"
  )
  expect_error(
    enrollment(lambda = Inf, N_total = 50),
    "finite positive"
  )
  expect_error(
    enrollment(lambda = 1, N_total = 0),
    "positive integer"
  )
  expect_error(
    enrollment(lambda = 1, N_total = 50.5),
    "positive integer"
  )
})

test_that("enrollment uses internal knots only", {
  expect_no_error(enrollment(lambda = 1, N_total = 5, lambda_time = NULL))
  expect_error(
    enrollment(lambda = c(1, 2), N_total = 50, lambda_time = NULL),
    "one greater"
  )
  expect_error(
    enrollment(lambda = 1, N_total = 50, lambda_time = 5),
    "one greater"
  )
  expect_error(
    enrollment(lambda = c(1, 2), N_total = 50, lambda_time = 0),
    "positive"
  )
  expect_error(
    enrollment(lambda = 1, N_total = 50, lambda_time = numeric()),
    "positive"
  )
})

test_that("enrollment validates all internal knots before simulation", {
  expect_error(
    enrollment(lambda = c(1, 2, 3), N_total = 50, lambda_time = c(20, 10)),
    "strictly increasing"
  )
  expect_error(
    enrollment(lambda = c(1, 2, 3), N_total = 50, lambda_time = c(10, 10)),
    "strictly increasing"
  )
  expect_error(
    enrollment(lambda = c(1, 2), N_total = 50, lambda_time = Inf),
    "finite positive"
  )
  expect_error(
    enrollment(lambda = c(1, 2), N_total = 50, lambda_time = -1),
    "finite positive"
  )
})
