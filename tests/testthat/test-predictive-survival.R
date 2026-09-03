predictive_survival_fixture <- function() {
  data.frame(
    time = c(2, 0, 4, 6, 0, 3),
    event = c(1L, 0L, 0L, 1L, 0L, 0L),
    treatment = c(0L, 1L, 1L, 0L, 0L, 1L),
    subject_enrolled = c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE),
    subject_impute_success = c(FALSE, FALSE, TRUE, FALSE, FALSE, TRUE),
    subject_impute_futility = c(FALSE, TRUE, FALSE, FALSE, TRUE, FALSE)
  )
}

predictive_survival_imputations <- function() {
  list(
    n_draws = 2L,
    current = list(
      rows = c(3L, 6L),
      time = matrix(c(7, 8, 5, 6), nrow = 2L),
      event = matrix(c(0L, 1L, 1L, 0L), nrow = 2L)
    ),
    future = list(
      rows = c(2L, 5L),
      time = matrix(c(4, 7, 8, 3), nrow = 2L),
      event = matrix(c(1L, 0L, 0L, 1L), nrow = 2L)
    )
  )
}

test_that("prepared predictive outcomes match completed patient data", {
  data <- predictive_survival_fixture()
  imputations <- predictive_survival_imputations()
  prepared <- goldilocks:::prepare_predictive_outcomes(
    data,
    imputations,
    check_futility = TRUE
  )

  for (draw in seq_len(imputations$n_draws)) {
    actual_current <- goldilocks:::complete_predictive_outcomes(
      prepared,
      imputations,
      draw
    )
    expected_current <- goldilocks:::complete_predictive_data(
      data,
      imputations,
      draw
    )
    expected_current <- expected_current[
      expected_current$subject_enrolled,
      c("time", "event", "treatment"),
      drop = FALSE
    ]

    expect_identical(actual_current$time, expected_current$time)
    expect_identical(actual_current$event, expected_current$event)
    expect_identical(actual_current$treatment, expected_current$treatment)

    actual_maximum <- goldilocks:::complete_predictive_outcomes(
      prepared,
      imputations,
      draw,
      maximum = TRUE
    )
    expected_maximum <- goldilocks:::complete_predictive_data(
      data,
      imputations,
      draw,
      include_future = TRUE
    )

    expect_identical(actual_maximum$time, expected_maximum$time)
    expect_identical(actual_maximum$event, expected_maximum$event)
    expect_identical(actual_maximum$treatment, expected_maximum$treatment)
  }
})

test_that("prepared outcomes omit an unused maximum sample", {
  data <- predictive_survival_fixture()
  imputations <- predictive_survival_imputations()
  prepared <- goldilocks:::prepare_predictive_outcomes(
    data,
    imputations,
    check_futility = FALSE
  )

  expect_null(prepared$maximum)
  expect_error(
    goldilocks:::complete_predictive_outcomes(
      prepared,
      imputations,
      draw = 1L,
      maximum = TRUE
    ),
    "maximum-sample outcomes unavailable"
  )
})

test_that("prepared sufficient statistics match the general route", {
  data <- data.frame(
    time = c(0, 4, 4.5, 8, 9, 12, 3, 7),
    event = c(0L, 1L, 0L, 1L, 0L, 1L, 1L, 0L),
    treatment = c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L)
  )
  rows <- c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE)

  checked <- goldilocks:::posterior_sufficient_stats(
    data = data,
    cutpoints = c(4, 8),
    single_arm = FALSE,
    rows = rows
  )
  prepared <- goldilocks:::summarise_survival_outcomes(
    time = data$time,
    event = data$event,
    treatment = data$treatment,
    cutpoints = c(4, 8),
    single_arm = FALSE,
    rows = rows
  )

  expect_identical(prepared, checked)
})

test_that("prepared survival analyses match completed data-frame analyses", {
  set.seed(3901)
  data <- data.frame(
    time = pmin(
      c(rexp(30, rate = 0.08), rexp(30, rate = 0.05)),
      12
    ),
    event = rep(c(1L, 1L, 0L), 20L),
    treatment = rep(0:1, each = 30L)
  )
  outcome <- list(
    time = data$time,
    event = data$event,
    treatment = data$treatment
  )
  normalized_prior <- goldilocks:::normalize_gamma_prior(
    prior_surv = c(0.1, 0.1),
    n_intervals = 3L,
    single_arm = FALSE,
    name = "prior_surv"
  )
  interval_widths <- goldilocks:::endpoint_interval_widths(c(4, 8), 12)

  for (method in c("logrank", "cox")) {
    expected <- goldilocks:::analyse_data(
      data = data,
      cutpoints = c(4, 8),
      end_of_study = 12,
      prior_surv = c(0.1, 0.1),
      N_mcmc = 1L,
      single_arm = FALSE,
      method = method,
      alternative = "less",
      h0 = 0
    )
    prepared <- goldilocks:::analyse_completed_survival(
      outcome = outcome,
      cutpoints = c(4, 8),
      end_of_study = 12,
      interval_widths = NULL,
      prior_surv = normalized_prior,
      N_mcmc = 1L,
      single_arm = FALSE,
      method = method,
      alternative = "less",
      h0 = 0,
      empty_interval = "prior"
    )

    expect_identical(prepared, expected, info = method)
  }

  set.seed(3902)
  expected_bayes <- goldilocks:::analyse_data(
    data = data,
    cutpoints = c(4, 8),
    end_of_study = 12,
    prior_surv = c(0.1, 0.1),
    N_mcmc = 100L,
    single_arm = FALSE,
    method = "bayes-surv",
    alternative = "less",
    h0 = 0,
    empty_interval = "prior"
  )
  expected_seed <- .Random.seed

  set.seed(3902)
  prepared_bayes <- goldilocks:::analyse_completed_survival(
    outcome = outcome,
    cutpoints = c(4, 8),
    end_of_study = 12,
    interval_widths = interval_widths,
    prior_surv = normalized_prior,
    N_mcmc = 100L,
    single_arm = FALSE,
    method = "bayes-surv",
    alternative = "less",
    h0 = 0,
    empty_interval = "prior"
  )

  expect_identical(prepared_bayes, expected_bayes)
  expect_identical(.Random.seed, expected_seed)
})
