test_that("pool_rubin_scalar applies Rubin's scalar pooling rules", {
  estimates <- c(-0.4, -0.2, -0.3)
  variances <- c(0.04, 0.05, 0.06)

  pooled <- pool_rubin_scalar(
    estimates = estimates,
    variances = variances,
    alternative = "two.sided",
    h0 = 0
  )

  expected_within <- mean(variances)
  expected_between <- var(estimates)
  expected_total <- expected_within + (1 + 1 / 3) * expected_between
  expected_relative <- (1 + 1 / 3) * expected_between / expected_within
  expected_df <- 2 * (1 + 1 / expected_relative)^2
  expected_statistic <- mean(estimates) / sqrt(expected_total)
  expected_success <- 1 - 2 * pt(-abs(expected_statistic), df = expected_df)

  expect_equal(pooled$estimate, mean(estimates))
  expect_equal(pooled$std_error, sqrt(expected_total))
  expect_equal(pooled$degrees_freedom, expected_df)
  expect_equal(pooled$success, expected_success)
})

test_that("pool_rubin_scalar preserves one-sided effect directions", {
  estimates <- c(-0.4, -0.2, -0.3)
  variances <- rep(0.05, 3)

  less <- pool_rubin_scalar(estimates, variances, "less", h0 = 0)
  greater <- pool_rubin_scalar(estimates, variances, "greater", h0 = 0)

  expect_gt(less$success, 0.5)
  expect_lt(greater$success, 0.5)
  expect_equal(less$success + greater$success, 1)
})

test_that("pool_rubin_scalar uses infinite degrees of freedom without between-imputation variation", {
  pooled <- pool_rubin_scalar(
    estimates = rep(-0.25, 3),
    variances = rep(0.04, 3),
    alternative = "less",
    h0 = 0
  )

  expect_equal(pooled$estimate, -0.25)
  expect_equal(pooled$std_error, 0.2)
  expect_identical(pooled$degrees_freedom, Inf)
  expect_equal(pooled$success, 1 - pnorm(-0.25 / 0.2))
})

test_that("pool_rubin_scalar can use between-imputation variance alone", {
  pooled <- pool_rubin_scalar(
    estimates = c(-0.1, 0, 0.1),
    variances = c(0, 0, 0),
    alternative = "two.sided",
    h0 = 0
  )

  expect_identical(pooled$degrees_freedom, 2)
  expect_gt(pooled$std_error, 0)
  expect_equal(pooled$estimate, 0)
})

test_that("prior_surv_final controls final-stage imputation", {
  data_in <- data.frame(
    time = rep(0.1, 20),
    event = 0,
    treatment = 1,
    subject_impute_success = TRUE
  )
  run_final <- function(prior_surv_final) {
    test_final(
      data_in = data_in,
      cutpoints = NULL,
      prior_surv_final = prior_surv_final,
      N_mcmc = 50,
      single_arm = TRUE,
      imputed_final = TRUE,
      method = "bayes-bin",
      N_impute = 5,
      alternative = "less",
      h0 = 0.5,
      prior_bin = c(1, 1),
      bin_method = "quadrature",
      binary_imputation = "event-time",
      empty_interval = "prior",
      end_of_study = 1
    )
  }

  set.seed(5109)
  low_hazard <- run_final(c(1e6, 1e14))
  set.seed(5109)
  high_hazard <- run_final(c(1e6, 1e4))

  expect_lt(low_hazard[2], 0.1)
  expect_gt(high_hazard[2], 0.9)
})
