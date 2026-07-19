test_that("prop_to_haz returns correct hazard for simple exponential", {
  lambda <- prop_to_haz(0.15, endtime = 36)
  expect_equal(pexp(36, lambda), 0.15)
})

test_that("prop_to_haz returns correct hazard for piecewise model", {
  haz <- prop_to_haz(c(0.15, 0.30), 12, 24)
  expect_length(haz, 2)
  # Verify using PWEALL
  p12 <- PWEALL::pwe(12, rate = haz, tchange = c(0, 12))$dist
  p24 <- PWEALL::pwe(24, rate = haz, tchange = c(0, 12))$dist
  expect_equal(p12, 0.15, tolerance = 1e-8)
  expect_equal(p24, 0.30, tolerance = 1e-8)
})

test_that("prop_to_haz returns positive hazard rates", {
  haz <- prop_to_haz(c(0.10, 0.25, 0.40), c(6, 12), 24)
  expect_true(all(haz > 0))
})

test_that("prop_to_haz errors on non-positive cutpoints", {
  expect_error(
    prop_to_haz(c(0.05, 0.15), cutpoints = 0, endtime = 36),
    "finite positive"
  )
})

test_that("prop_to_haz validates cumulative probabilities", {
  expect_error(
    prop_to_haz(c(0.30, 0.20), 12, 24),
    "non-decreasing"
  )

  expect_error(
    prop_to_haz(1, cutpoints = NULL, endtime = 36),
    "probs"
  )

  expect_error(
    prop_to_haz(1.1, cutpoints = NULL, endtime = 36),
    "probs"
  )

  expect_error(
    prop_to_haz(NA_real_, cutpoints = NULL, endtime = 36),
    "probs"
  )
})

test_that("prop_to_haz validates cutpoint and endpoint dimensions", {
  expect_error(
    prop_to_haz(c(0.15, 0.30), cutpoints = NULL, endtime = 36),
    "one greater"
  )

  expect_error(
    prop_to_haz(c(0.15, 0.30), cutpoints = 12, endtime = 12),
    "greater than the last cutpoint"
  )

  expect_error(
    prop_to_haz(c(0.10, 0.20, 0.30), cutpoints = c(12, 12), endtime = 36),
    "strictly increasing"
  )
})
