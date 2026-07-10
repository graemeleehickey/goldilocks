test_that("enrollment returns correct length", {
  set.seed(7413)
  out <- enrollment(lambda = 5, N_total = 50, lambda_time = 0)
  expect_length(out, 50)
})

test_that("enrollment times are non-decreasing and start at zero", {
  set.seed(2198)
  out <- enrollment(lambda = 5, N_total = 100, lambda_time = 0)
  expect_equal(out[1], 0)
  expect_false(is.unsorted(out))
})

test_that("enrollment works with piecewise rates", {
  set.seed(4821)
  out <- enrollment(
    lambda = c(0.3, 0.7, 1.2),
    N_total = 80,
    lambda_time = c(0, 10, 20)
  )
  expect_length(out, 80)
  expect_equal(out[1], 0)
  expect_false(is.unsorted(out))
})

test_that("enrollment errors on non-positive lambda", {
  expect_error(
    enrollment(lambda = -1, N_total = 50, lambda_time = 0),
    "finite positive"
  )
  expect_error(
    enrollment(lambda = 0, N_total = 50, lambda_time = 0),
    "finite positive"
  )
})

test_that("enrollment errors on non-positive N_total", {
  expect_error(
    enrollment(lambda = 1, N_total = 0, lambda_time = 0),
    "positive integer"
  )
})

test_that("enrollment errors on mismatched lambda and lambda_time lengths", {
  expect_error(
    enrollment(lambda = c(1, 2), N_total = 50, lambda_time = 0),
    "length"
  )
})

test_that("enrollment errors when lambda_time does not start at 0", {
  expect_error(
    enrollment(lambda = 1, N_total = 50, lambda_time = 5),
    "first cutpoint"
  )
})

test_that("enrollment uses correct rate at changepoint boundaries", {
  # Use very different rates so the boundary behavior is detectable:
  # rate 1 before time 50, rate 100 after time 50.
  # With the old >= bug, time 50 would still use rate 1.
  set.seed(6739)
  out <- enrollment(
    lambda = c(1, 100),
    N_total = 500,
    lambda_time = c(0, 50)
  )
  # Most subjects should enroll at or after the changepoint since the
  # rate jumps 100x. Check that the median enrollment time is near 50.
  expect_true(median(out) >= 45)
})

test_that("enrollment errors on NULL lambda_time", {
  expect_error(
    enrollment(lambda = 1, N_total = 50, lambda_time = NULL),
    "length"
  )
})

test_that("enrollment validates all schedule inputs before simulation", {
  expect_error(
    enrollment(lambda = Inf, N_total = 50, lambda_time = 0),
    "finite positive"
  )
  expect_error(
    enrollment(lambda = 1, N_total = 50.5, lambda_time = 0),
    "positive integer"
  )
  expect_error(
    enrollment(lambda = c(1, 2, 3), N_total = 50, lambda_time = c(0, 20, 10)),
    "strictly increasing"
  )
  expect_error(
    enrollment(lambda = c(1, 2), N_total = 50, lambda_time = c(0, 0)),
    "strictly increasing"
  )
  expect_error(
    enrollment(lambda = c(1, 2), N_total = 50, lambda_time = c(0, Inf)),
    "finite numeric"
  )
})
