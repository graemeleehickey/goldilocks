test_that("prop_to_haz returns correct hazard for simple exponential", {
  lambda <- prop_to_haz(0.15, endtime = 36)
  expect_equal(pexp(36, lambda), 0.15)
})

test_that("prop_to_haz returns correct hazard for piecewise model", {
  haz <- prop_to_haz(c(0.15, 0.30), c(0, 12), 24)
  expect_length(haz, 2)
  # Verify using PWEALL
  p12 <- PWEALL::pwe(12, rate = haz, tchange = c(0, 12))$dist
  p24 <- PWEALL::pwe(24, rate = haz, tchange = c(0, 12))$dist
  expect_equal(p12, 0.15, tolerance = 1e-8)
  expect_equal(p24, 0.30, tolerance = 1e-8)
})

test_that("prop_to_haz returns positive hazard rates", {
  haz <- prop_to_haz(c(0.10, 0.25, 0.40), c(0, 6, 12), 24)
  expect_true(all(haz > 0))
})

test_that("prop_to_haz errors when cutpoints don't start at 0", {
  expect_error(
    prop_to_haz(0.15, cutpoints = 5, endtime = 36),
    "First element"
  )
})
