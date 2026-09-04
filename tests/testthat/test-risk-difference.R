test_that("Farrington-Manning statistics agree with reference values", {
  # Independently generated with gsDesign::testBinomial() 3.9.0 using
  # scale = "Difference" and no Miettinen-Nurminen adjustment.
  cases <- data.frame(
    events_treatment = c(20, 3, 9, 56),
    events_control = c(10, 1, 4, 48),
    n_treatment = c(50, 40, 30, 70),
    n_control = c(50, 25, 20, 80),
    h0 = c(0, 0.05, -0.1, 0.1),
    expected = c(
      2.18217890235992,
      -0.263824635862046,
      1.53525281484015,
      1.34773871811343
    )
  )

  observed <- vapply(
    seq_len(nrow(cases)),
    function(index) {
      risk_difference_fm_from_counts(
        events_control = cases$events_control[[index]],
        n_control = cases$n_control[[index]],
        events_treatment = cases$events_treatment[[index]],
        n_treatment = cases$n_treatment[[index]],
        alternative = "greater",
        h0 = cases$h0[[index]]
      )$statistic
    },
    numeric(1)
  )

  expect_equal(observed, cases$expected, tolerance = 1e-7)
})

test_that("Farrington-Manning handles sparse boundary tables", {
  all_zero <- risk_difference_fm_from_counts(
    events_control = 0,
    n_control = 20,
    events_treatment = 0,
    n_treatment = 20,
    alternative = "greater",
    h0 = 0
  )
  all_one <- risk_difference_fm_from_counts(
    events_control = 20,
    n_control = 20,
    events_treatment = 20,
    n_treatment = 20,
    alternative = "greater",
    h0 = 0
  )
  treatment_only <- risk_difference_fm_from_counts(
    events_control = 0,
    n_control = 20,
    events_treatment = 1,
    n_treatment = 20,
    alternative = "greater",
    h0 = 0
  )

  for (fit in list(all_zero, all_one)) {
    expect_identical(fit$statistic, 0)
    expect_identical(fit$std_error, 0)
    expect_equal(fit$success, 0.5)
  }
  expect_true(is.finite(treatment_only$statistic))
  expect_gt(treatment_only$success, 0.5)

  all_zero_two_sided <- risk_difference_fm_from_counts(
    events_control = 0,
    n_control = 20,
    events_treatment = 0,
    n_treatment = 20,
    alternative = "two.sided",
    h0 = 0
  )
  expect_identical(all_zero_two_sided$success, 0)
})

test_that("Farrington-Manning preserves direction and arm-swap symmetry", {
  greater <- risk_difference_fm_from_counts(
    events_control = 2,
    n_control = 30,
    events_treatment = 9,
    n_treatment = 40,
    alternative = "greater",
    h0 = 0.1
  )
  less <- risk_difference_fm_from_counts(
    events_control = 2,
    n_control = 30,
    events_treatment = 9,
    n_treatment = 40,
    alternative = "less",
    h0 = 0.1
  )
  swapped <- risk_difference_fm_from_counts(
    events_control = 9,
    n_control = 40,
    events_treatment = 2,
    n_treatment = 30,
    alternative = "less",
    h0 = -0.1
  )

  expect_equal(greater$success + less$success, 1)
  expect_equal(swapped$statistic, -greater$statistic)
  expect_equal(swapped$success, greater$success)
  expect_equal(swapped$estimate, -greater$estimate)
})

test_that("Farrington-Manning has stable sparse-null rejection rates", {
  exact_rejection_rate <- function(
    control_risk,
    treatment_risk,
    n_control,
    n_treatment,
    h0
  ) {
    outcomes <- expand.grid(
      events_control = 0:n_control,
      events_treatment = 0:n_treatment
    )
    success <- mapply(
      function(events_control, events_treatment) {
        risk_difference_fm_from_counts(
          events_control = events_control,
          n_control = n_control,
          events_treatment = events_treatment,
          n_treatment = n_treatment,
          alternative = "greater",
          h0 = h0
        )$success
      },
      outcomes$events_control,
      outcomes$events_treatment
    )
    probability <-
      dbinom(outcomes$events_control, n_control, control_risk) *
      dbinom(outcomes$events_treatment, n_treatment, treatment_risk)
    sum(probability[success > 0.975])
  }

  rates <- c(
    equal_rare_20 = exact_rejection_rate(0.05, 0.05, 20, 20, 0),
    equal_rare_50 = exact_rejection_rate(0.01, 0.01, 50, 50, 0),
    nonzero_null = exact_rejection_rate(0.15, 0.05, 40, 60, -0.1)
  )

  expect_equal(
    rates,
    c(
      equal_rare_20 = 0.0058252846,
      equal_rare_50 = 0.0009659066,
      nonzero_null = 0.0194879064
    ),
    tolerance = 1e-9
  )
  expect_true(all(rates <= 0.025))
})

test_that("riskdiff-fm is available through the completed-data analysis", {
  data <- data.frame(
    time = rep(36, 40),
    event = c(rep(0, 20), 1, rep(0, 19)),
    treatment = rep(0:1, each = 20)
  )

  result <- analyse_data(
    data = data,
    cutpoints = NULL,
    end_of_study = 36,
    prior_surv = c(0.1, 0.1),
    N_mcmc = 1,
    single_arm = FALSE,
    method = "riskdiff-fm",
    alternative = "greater",
    h0 = 0
  )

  expect_named(result, c("success", "effect"))
  expect_gt(result$success, 0.5)
  expect_equal(result$effect, 0.05)
})

test_that("explicit Wald analysis retains its zero-variance error", {
  fit <- risk_difference_from_counts(
    events_control = 0,
    n_control = 20,
    events_treatment = 0,
    n_treatment = 20
  )

  expect_error(
    risk_difference_wald_from_estimate(fit, "greater", h0 = 0),
    "estimated variance is zero"
  )
})

test_that("sparse nonzero-null simulations use riskdiff-fm end to end", {
  out <- sim_trials(
    hazard_treatment = prop_to_haz(0.10, endtime = 1),
    hazard_control = prop_to_haz(0.05, endtime = 1),
    N_total = 40,
    lambda = 100,
    end_of_study = 1,
    alternative = "greater",
    h0 = 0.05,
    prob_ha = 0.975,
    N_impute = 1,
    N_mcmc = 1,
    N_trials = 6,
    method = "riskdiff-fm",
    backend = "sequential",
    seed = 6363
  )

  expect_equal(
    out$sims$post_prob_ha,
    c(0.5, 0.15245089, 0.5, 0.83946882, 0.5, 0.30261104),
    tolerance = 1e-8
  )
  expect_equal(out$sims$est_final, c(0.05, 0, 0.05, 0.15, 0.05, 0))
  expect_equal(nrow(out$failures), 0)
  expect_identical(attr(out, "arguments")$method, "riskdiff-fm")
})

test_that("sim_trials warns once for the deprecated riskdiff alias", {
  warnings <- character()
  out <- withCallingHandlers(
    sim_trials(
      hazard_treatment = prop_to_haz(0.2, endtime = 1),
      hazard_control = prop_to_haz(0.4, endtime = 1),
      N_total = 40,
      lambda = 100,
      end_of_study = 1,
      alternative = "less",
      N_trials = 2,
      method = "riskdiff",
      backend = "sequential",
      seed = 6364
    ),
    warning = function(warning) {
      warnings <<- c(warnings, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  )

  expect_length(warnings, 1)
  expect_match(warnings, "deprecated; use `method = \"riskdiff-wald\"`")
  expect_identical(attr(out, "arguments")$method, "riskdiff-wald")
})
