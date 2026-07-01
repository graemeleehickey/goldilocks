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
})
