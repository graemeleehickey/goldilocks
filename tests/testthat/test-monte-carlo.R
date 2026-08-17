test_that("finite Monte Carlo bounds control expected-success decision error", {
  draws <- 200L
  threshold <- 0.9
  summaries <- lapply(0:draws, function(successes) {
    monte_carlo_probability_summary(
      successes,
      draws,
      threshold,
      direction = "greater",
      confidence = 0.95
    )
  })
  crossing_counts <- which(vapply(
    summaries,
    function(x) x$crossed,
    logical(1)
  )) - 1L
  decision_error <- sum(dbinom(crossing_counts, draws, threshold))

  expect_lte(decision_error, 0.05)
})

test_that("finite Monte Carlo bounds control futility decision error", {
  draws <- 200L
  threshold <- 0.05
  summaries <- lapply(0:draws, function(successes) {
    monte_carlo_probability_summary(
      successes,
      draws,
      threshold,
      direction = "less",
      confidence = 0.95
    )
  })
  crossing_counts <- which(vapply(
    summaries,
    function(x) x$crossed,
    logical(1)
  )) - 1L
  decision_error <- sum(dbinom(crossing_counts, draws, threshold))

  expect_lte(decision_error, 0.05)
})

test_that("coarse and equality-bound Monte Carlo results continue", {
  coarse <- monte_carlo_probability_summary(
    successes = 10,
    draws = 10,
    threshold = 0.9,
    direction = "greater"
  )
  equality <- monte_carlo_probability_summary(
    successes = 90,
    draws = 100,
    threshold = 0.9,
    direction = "greater"
  )

  expect_true(coarse$point_crossed)
  expect_false(coarse$crossed)
  expect_identical(coarse$reason, "monte_carlo_uncertain")
  expect_false(equality$point_crossed)
  expect_false(equality$crossed)
})

test_that("inner posterior Monte Carlo uncertainty is not called success", {
  analysis <- set_analysis_mc_counts(
    list(success = 1, effect = 0.2),
    successes = 10,
    draws = 10
  )
  classification <- classify_completed_analysis(
    analysis,
    prob_ha = 0.95,
    mc_conf_level = 0.95
  )

  expect_false(classification$crossed)
  expect_true(classification$uncertain)
})

test_that("Bayesian completed-data analyses retain their Monte Carlo counts", {
  data <- data.frame(
    time = rep(12, 20),
    event = c(rep(0, 10), rep(1, 10)),
    treatment = rep(0:1, each = 10)
  )
  result <- analyse_data(
    data = data,
    cutpoints = NULL,
    end_of_study = 12,
    prior_surv = c(0.1, 0.1),
    N_mcmc = 25,
    single_arm = FALSE,
    method = "bayes-bin",
    alternative = "greater",
    h0 = 0,
    prior_bin = c(1, 1),
    bin_method = "mc"
  )
  counts <- attr(result, "mc_counts", exact = TRUE)

  expect_identical(counts$draws, 25L)
  expect_gte(counts$successes, 0L)
  expect_lte(counts$successes, counts$draws)
})

test_that("Monte Carlo defaults no longer imply ten-draw decisions", {
  expect_identical(formals(survival_adapt)$N_impute, 500)
  expect_identical(formals(survival_adapt)$N_mcmc, 1000)
  expect_identical(formals(sim_trials)$N_impute, 500)
  expect_identical(formals(sim_trials)$N_mcmc, 1000)
})
