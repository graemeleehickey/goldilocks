test_that("posterior returns correct dimensions (two-arm, single interval)", {
  set.seed(3821)
  data <- data.frame(
    time      = rexp(80, 0.02),
    event     = sample(0:1, 80, replace = TRUE),
    treatment = rep(0:1, each = 40)
  )
  res <- posterior(data, cutpoints = 0, prior = c(0.1, 0.1),
                   N_mcmc = 50, single_arm = FALSE)
  expect_true(is.array(res))
  expect_equal(dim(res), c(50, 1, 2))
  expect_true(all(res[, , 1] > 0))
  expect_true(all(res[, , 2] > 0))
})

test_that("posterior returns correct dimensions (two-arm, piecewise)", {
  set.seed(6194)
  data <- data.frame(
    time      = rexp(100, 0.02),
    event     = sample(0:1, 100, replace = TRUE),
    treatment = rep(0:1, each = 50)
  )
  res <- posterior(data, cutpoints = c(0, 12, 24), prior = c(0.1, 0.1),
                   N_mcmc = 30, single_arm = FALSE)
  expect_equal(dim(res), c(30, 3, 2))
})

test_that("posterior returns correct dimensions (single-arm)", {
  set.seed(7452)
  data <- data.frame(
    time      = rexp(50, 0.02),
    event     = sample(0:1, 50, replace = TRUE),
    treatment = rep(1, 50)
  )
  res <- posterior(data, cutpoints = 0, prior = c(0.1, 0.1),
                   N_mcmc = 40, single_arm = TRUE)
  expect_equal(dim(res), c(40, 1, 2))
  expect_true(all(res[, , 1] > 0))
  # Control slice should be all NA for single-arm
  expect_true(all(is.na(res[, , 2])))
})

test_that("posterior warns on zero-exposure interval", {
  set.seed(2913)
  # All event times < 10, so interval 2 (>= 10) has zero exposure
  data <- data.frame(
    time      = runif(60, 0, 9),
    event     = sample(0:1, 60, replace = TRUE),
    treatment = rep(0:1, each = 30)
  )
  expect_warning(
    posterior(data, cutpoints = c(0, 10), prior = c(0.1, 0.1),
             N_mcmc = 10, single_arm = FALSE),
    "zero subjects"
  )
})

test_that("posterior returns positive draws", {
  set.seed(5087)
  data <- data.frame(
    time      = rexp(100, 0.03),
    event     = sample(0:1, 100, replace = TRUE),
    treatment = rep(0:1, each = 50)
  )
  res <- posterior(data, cutpoints = c(0, 12), prior = c(0.1, 0.1),
                   N_mcmc = 100, single_arm = FALSE)
  expect_true(all(res > 0))
})
