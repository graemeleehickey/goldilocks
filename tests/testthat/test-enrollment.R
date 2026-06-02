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
    "non-negative"
  )
  expect_error(
    enrollment(lambda = 0, N_total = 50, lambda_time = 0),
    "non-negative"
  )
})

test_that("enrollment errors on non-positive N_total", {
  expect_error(
    enrollment(lambda = 1, N_total = 0, lambda_time = 0),
    "greater than 0"
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

test_that("enrollment errors on NULL lambda_time", {
  expect_error(
    enrollment(lambda = 1, N_total = 50, lambda_time = NULL),
    "length"
  )
})
