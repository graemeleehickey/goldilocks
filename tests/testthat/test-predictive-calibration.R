near_futility_interim_fixture <- function() {
  data.frame(
    id = 1:8,
    treatment = c(0, 1, 0, 1, 0, 1, 0, 1),
    enrollment = 0:7,
    time = c(8, 7, 6, 5, 4, 3, 2, 1),
    event = c(1, 0, 0, 1, 0, 0, 0, 0),
    status = c(
      "event",
      "pending",
      "censored",
      "event",
      "pending",
      "pending",
      "pending",
      "pending"
    ),
    stringsAsFactors = FALSE
  )
}

near_success_interim_fixture <- function() {
  n <- 20L
  treatment <- rep(0:1, length.out = n)
  enrollment <- 0:(n - 1L)
  available_followup <- pmin(12, 20 - enrollment)
  control_events <- head(which(treatment == 0), 5L)
  treatment_event <- head(which(treatment == 1), 1L)
  event_ids <- c(control_events, treatment_event)
  event <- as.integer(seq_len(n) %in% event_ids)
  time <- available_followup
  time[control_events] <- 0.55 * available_followup[control_events]
  time[treatment_event] <- 6.7
  status <- ifelse(
    event == 1,
    "event",
    ifelse(available_followup >= 12, "complete", "pending")
  )

  data.frame(
    id = seq_len(n),
    treatment = treatment,
    enrollment = enrollment,
    time = time,
    event = event,
    status = status,
    stringsAsFactors = FALSE
  )
}

test_that("maximum-sample prediction is calibrated near the futility threshold", {
  skip_on_cran()

  # Independent reference: seed 9801, 250,000 imputations.
  reference_successes <- 12787L
  reference_draws <- 250000L
  reference_estimate <- reference_successes / reference_draws
  reference_mcse <- sqrt(
    reference_estimate * (1 - reference_estimate) / reference_draws
  )

  out <- evaluate_interim(
    data = near_futility_interim_fixture(),
    data_cut = 8,
    look = 2,
    N_total = 12,
    end_of_study = 10,
    rand_ratio = c(control = 1, treatment = 1),
    alternative = "less",
    Fn = 0.05,
    Sn = 0.9,
    prob_ha = 0.95,
    N_impute = 50000,
    N_mcmc = 1,
    method = "logrank",
    seed = 9811
  )

  estimate <- out$monte_carlo[
    out$monte_carlo$estimand == "success_at_maximum",
    ,
    drop = FALSE
  ]
  probability <- out$probabilities[
    out$probabilities$estimand == "success_at_maximum",
    ,
    drop = FALSE
  ]

  expect_lt(abs(reference_estimate - 0.05), 0.005)
  expect_mc_close(
    estimate$estimate,
    reference_estimate,
    estimate$mcse,
    "predictive success at maximum sample size near futility",
    reference_mcse = reference_mcse
  )
  expect_identical(
    probability$threshold_crossed,
    probability$probability < probability$threshold
  )
})

test_that("current-sample prediction is calibrated near the success threshold", {
  skip_on_cran()

  # Independent reference: seed 9802, 250,000 imputations.
  reference_successes <- 225621L
  reference_draws <- 250000L
  reference_estimate <- reference_successes / reference_draws
  reference_mcse <- sqrt(
    reference_estimate * (1 - reference_estimate) / reference_draws
  )

  out <- evaluate_interim(
    data = near_success_interim_fixture(),
    data_cut = 20,
    look = 1,
    N_total = 30,
    end_of_study = 12,
    rand_ratio = c(control = 1, treatment = 1),
    alternative = "less",
    Fn = 0.05,
    Sn = 0.9,
    prob_ha = 0.95,
    N_impute = 50000,
    N_mcmc = 1,
    method = "logrank",
    seed = 9812
  )

  estimate <- out$monte_carlo[
    out$monte_carlo$estimand == "success_if_stop_now",
    ,
    drop = FALSE
  ]
  probability <- out$probabilities[
    out$probabilities$estimand == "success_if_stop_now",
    ,
    drop = FALSE
  ]

  expect_lt(abs(reference_estimate - 0.9), 0.005)
  expect_mc_close(
    estimate$estimate,
    reference_estimate,
    estimate$mcse,
    "predictive success if accrual stops near expected-success threshold",
    reference_mcse = reference_mcse
  )
  expect_identical(
    probability$threshold_crossed,
    probability$probability > probability$threshold
  )
})
