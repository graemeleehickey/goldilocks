test_that("dropout censors only when it precedes the event and horizon", {
  local_mocked_bindings(
    enrollment = function(...) 0:5,
    pwe_sim = function(...) {
      data.frame(time = c(2, 8, 10, 10, 5, 10), event = c(1, 1, 0, 0, 1, 0))
    },
    rexp = function(n, rate) {
      expect_equal(n, 6)
      expect_equal(rate, -log(0.8) / 10)
      c(1, 9, 4, 20, 5, 10)
    }
  )
  out <- sim_comp_data(
    hazard_treatment = 0.1,
    N_total = 6,
    end_of_study = 10,
    prop_loss = 0.2
  )
  expect_equal(out$time, c(1, 8, 4, 10, 5, 10))
  expect_equal(out$event, c(0, 1, 0, 0, 1, 0))
  expect_identical(out$loss_to_fu, c(TRUE, FALSE, TRUE, FALSE, FALSE, FALSE))
})

test_that("zero dropout preserves event data and consumes no extra RNG", {
  set.seed(6101)
  expected_enrollment <- enrollment(lambda = 5, N_total = 40)
  expected_treatment <- randomization(N_total = 40, block = 2)
  control <- pwe_sim(n = 20, hazard = 0.02, maxtime = 36)
  treatment <- pwe_sim(n = 20, hazard = 0.01, maxtime = 36)
  expected_time <- expected_event <- numeric(40)
  expected_time[expected_treatment == 0] <- control$time
  expected_time[expected_treatment == 1] <- treatment$time
  expected_event[expected_treatment == 0] <- control$event
  expected_event[expected_treatment == 1] <- treatment$event
  expected_seed <- .Random.seed

  for (p in list(0, c(treatment = 0, control = 0))) {
    set.seed(6101)
    out <- sim_comp_data(
      hazard_treatment = 0.01,
      hazard_control = 0.02,
      N_total = 40,
      lambda = 5,
      end_of_study = 36,
      prop_loss = p
    )
    expect_identical(.Random.seed, expected_seed)
    expect_equal(out$enrollment, expected_enrollment)
    expect_equal(out$treatment, expected_treatment)
    expect_identical(out$time, expected_time)
    expect_identical(out$event, expected_event)
    expect_false(any(out$loss_to_fu))
  }
})

test_that("dropout CDF is calibrated to the per-subject follow-up horizon", {
  # With no endpoint events, actual dropout directly reveals the specified
  # dropout-time CDF. Check an earlier time as well as the reference horizon.
  set.seed(6102)
  out <- sim_comp_data(
    hazard_treatment = 0,
    hazard_control = 0,
    N_total = 40000,
    lambda = 5,
    end_of_study = 12,
    prop_loss = c(control = 0, treatment = 0.3)
  )
  control <- out[out$treatment == 0, ]
  treatment <- out[out$treatment == 1, ]
  expect_false(any(control$loss_to_fu))
  expect_true(all(control$time == 12))
  for (t in c(6, 12)) {
    expect_mc_rate(
      treatment$loss_to_fu & treatment$time <= t,
      target = 1 - 0.7^(t / 12),
      estimand = paste("dropout CDF at month", t)
    )
  }
  expect_true(all(out$event == 0))
})

test_that("small trials have random loss counts without forced rounding", {
  set.seed(6103)
  # Previously every one-patient trial was forced to lose its only patient.
  lost <- replicate(2000, {
    sim_comp_data(
      hazard_treatment = 0,
      N_total = 1,
      end_of_study = 12,
      prop_loss = 0.05
    )$loss_to_fu
  })
  expect_mc_rate(lost, 0.05, "dropout probability in one-patient trials")
})

test_that("independent dropout recovers survival and observed outcome probabilities", {
  scenarios <- list(
    single = list(
      hazard_treatment = log(2) / 36,
      prop_loss = 0.3
    ),
    common = list(
      hazard_treatment = 0.01,
      hazard_control = 0.02,
      prop_loss = 0.3
    ),
    differential_piecewise = list(
      hazard_treatment = c(0.03, 0.01),
      hazard_control = c(0.02, 0.04),
      generation_cutpoints = 12,
      prop_loss = c(control = 0.1, treatment = 0.4)
    )
  )
  horizon <- 36
  for (scenario in names(scenarios)) {
    args <- scenarios[[scenario]]
    set.seed(6104)
    out <- do.call(
      sim_comp_data,
      c(
        args,
        list(
          N_total = 40000,
          lambda = 20,
          end_of_study = horizon
        )
      )
    )
    dropout <- normalize_prop_loss(args$prop_loss, is.null(args$hazard_control))
    for (arm in names(dropout)) {
      data <- out[out$treatment == if (arm == "treatment") 1 else 0, ]
      hazard <- args[[paste0("hazard_", arm)]]
      starts <- c(0, args$generation_cutpoints)
      widths <- diff(c(starts, horizon))
      eta <- -log1p(-dropout[[arm]]) / horizon
      total_rate <- hazard + eta
      # Integrate the two independent cause-specific densities over each
      # generating interval, including survival through earlier intervals.
      survival_start <- exp(-c(0, head(cumsum(total_rate * widths), -1)))
      exposure <- survival_start * (-expm1(-total_rate * widths)) / total_rate
      label <- paste(scenario, arm)
      expect_mc_rate(
        data$loss_to_fu,
        sum(eta * exposure),
        paste(label, "observed dropout")
      )
      expect_mc_rate(
        data$event == 1,
        sum(hazard * exposure),
        paste(label, "observed events")
      )

      times <- c(12, 24, horizon)
      fit <- summary(
        survival::survfit(survival::Surv(time, event) ~ 1, data = data),
        times = times
      )
      target <- vapply(
        times,
        function(t) {
          exp(-sum(hazard * pmax(0, pmin(t - starts, widths))))
        },
        numeric(1)
      )
      for (j in seq_along(times)) {
        expect_mc_close(
          estimate = fit$surv[j],
          target = target[j],
          mcse = fit$std.err[j],
          estimand = paste(label, "Kaplan-Meier survival at", times[j])
        )
      }
    }
  }
})

test_that("dropout probabilities exclude one and reject malformed inputs", {
  invalid <- list(
    1,
    c(control = 0.1, treatment = 1),
    -0.1,
    1.1,
    NA_real_,
    NaN,
    Inf,
    numeric(),
    "0.1",
    matrix(0.1)
  )
  for (p in invalid) {
    expect_error(
      sim_comp_data(
        hazard_treatment = 0.01,
        hazard_control = 0.02,
        N_total = 20,
        end_of_study = 12,
        prop_loss = p
      ),
      "prop_loss.*finite probabilities.*\\[0, 1\\)"
    )
  }
  # Validate before expensive trial simulation, including both public wrappers.
  expect_error(
    survival_adapt(
      hazard_treatment = 0.01,
      N_total = 20,
      end_of_study = 12,
      prop_loss = 1,
      hazard_control = NULL
    ),
    "prop_loss.*\\[0, 1\\)"
  )
  expect_error(
    sim_trials(
      hazard_treatment = 0.01,
      N_total = 20,
      end_of_study = 12,
      prop_loss = 1,
      hazard_control = NULL,
      N_trials = 1
    ),
    "prop_loss.*\\[0, 1\\)"
  )
  set.seed(6105)
  for (p in c(1e-12, 1 - .Machine$double.eps)) {
    out <- sim_comp_data(
      hazard_treatment = 0.01,
      N_total = 20,
      end_of_study = 12,
      prop_loss = p
    )
    expect_true(all(is.finite(out$time) & out$time > 0 & out$time <= 12))
  }
})
