test_that("calendar metrics use first enrollment as their time origin", {
  data <- data.frame(
    enrollment = c(0, 2, 4),
    time = c(5, 2, 1)
  )

  metrics <- trial_calendar_metrics(data, end_of_study = 6)

  expect_equal(metrics$accrual_stop_time, 4)
  expect_equal(metrics$analysis_ready_time, 5)
  expect_equal(metrics$planned_completion_time, 10)
  expect_equal(metrics$followup_person_time, 8)
  expect_equal(metrics$peak_active_followup, 2L)
  expect_equal(active_followup_at(data, 2), 2)
  expect_equal(active_followup_at(data, 4), 2)
  expect_equal(active_followup_at(data, 5), 0)
})

test_that("analysis readiness cannot precede the end of accrual", {
  data <- data.frame(
    enrollment = c(0, 10),
    time = c(1, 0)
  )

  metrics <- trial_calendar_metrics(data, end_of_study = 12)

  expect_equal(metrics$analysis_ready_time, 10)
  expect_equal(metrics$peak_active_followup, 1L)
})

calendar_test_result <- function(include_traces = TRUE) {
  sims <- data.frame(
    N_enrolled = c(100, 150, 200),
    stop_futility = c(0, 1, 0),
    stop_expected_success = c(1, 0, 0),
    stopping_reason = c(
      "expected_success",
      "futility",
      "maximum_sample_size"
    ),
    accrual_stop_time = c(10, 15, 20),
    analysis_ready_time = c(20, 25, 30),
    planned_completion_time = c(22, 27, 32),
    followup_person_time = c(800, 1100, 1500),
    peak_active_followup = c(50, 60, 70)
  )
  failures <- data.frame(
    trial = 4L,
    error_class = "simpleError",
    message = "synthetic failure"
  )
  out <- list(
    sims = sims,
    failures = failures,
    call = quote(sim_trials(N_trials = 4))
  )
  if (include_traces) {
    out$traces <- data.frame(
      trial = c(1L, 2L, 3L, 2L, 3L),
      look = c(1L, 1L, 1L, 2L, 2L),
      planned_N = c(100L, 100L, 100L, 150L, 150L),
      calendar_time = c(9, 10, 11, 14, 16),
      active_followup = c(40L, 42L, 44L, 50L, 54L)
    )
  }
  attr(out, "parallel_metadata") <- list(
    backend = "sequential",
    tasks = 4L
  )
  out
}

test_that("calendar summaries are wide and retain failed denominators", {
  out <- summarise_calendar_time(calendar_test_result())

  expect_s3_class(out, "goldilocks_calendar_summary")
  expect_named(out, c("trial_duration", "interim_timing"))

  duration <- out$trial_duration
  expect_equal(
    duration$stopping_reason,
    c("expected_success", "futility", "maximum_sample_size", "overall")
  )
  expect_equal(duration$n_requested, rep(4L, 4))
  expect_equal(duration$n_analyzed, rep(3L, 4))
  expect_equal(duration$n_failed, rep(1L, 4))
  expect_equal(duration$n_trials, c(1L, 1L, 1L, 3L))
  expect_equal(duration$percent_trials, c(25, 25, 25, 75))
  expect_equal(duration$analysis_ready_median, c(20, 25, 30, 25))
  expect_equal(duration$followup_person_time_mean[4], 3400 / 3)
  expect_true(all(
    c(
      "accrual_stop_se",
      "analysis_ready_p10",
      "planned_completion_p90",
      "peak_active_followup_mean"
    ) %in%
      names(duration)
  ))
  expect_false(any(grepl("^decision_time_", names(duration))))

  interim <- out$interim_timing
  expect_equal(interim$look, c(1L, 2L))
  expect_equal(interim$n_requested, c(4L, 4L))
  expect_equal(interim$n_analyzed, c(3L, 3L))
  expect_equal(interim$n_failed, c(1L, 1L))
  expect_equal(interim$n_reached, c(3L, 2L))
  expect_equal(interim$percent_reached, c(75, 50))
  expect_equal(interim$calendar_time_median, c(10, 15))
  expect_equal(interim$active_followup_mean, c(42, 52))
})

test_that("calendar summaries distinguish immediate decision timing", {
  sims <- data.frame(
    N_enrolled = c(100, 120, 140, 200),
    stop_futility = c(FALSE, FALSE, TRUE, FALSE),
    stop_immediate_success = c(TRUE, FALSE, FALSE, FALSE),
    stop_expected_success = c(FALSE, TRUE, FALSE, FALSE),
    accrual_stop_time = c(10, 12, 14, 20),
    decision_time = c(10, 22, 14, 30),
    analysis_ready_time = c(20, 22, 24, 30),
    planned_completion_time = c(22, 24, 26, 32),
    followup_person_time = c(800, 900, 1000, 1500),
    peak_active_followup = c(50, 55, 60, 70)
  )

  expect_warning(
    out <- summarise_calendar_time(sims),
    "return_trace = TRUE"
  )
  duration <- out$trial_duration

  expect_equal(
    duration$stopping_reason,
    c(
      "immediate_success",
      "expected_success",
      "futility",
      "maximum_sample_size",
      "overall"
    )
  )
  expect_equal(duration$decision_time_median, c(10, 22, 14, 30, 18))
  expect_true(all(
    c(
      "decision_time_mean",
      "decision_time_se",
      "decision_time_p10",
      "decision_time_p90"
    ) %in%
      names(duration)
  ))
})

test_that("calendar stopping reason gives immediate success precedence", {
  expect_equal(
    calendar_stopping_reason(
      stop_futility = c(TRUE, FALSE, FALSE, FALSE),
      stop_expected_success = c(FALSE, TRUE, FALSE, FALSE),
      stop_immediate_success = c(TRUE, FALSE, FALSE, FALSE)
    ),
    c(
      "immediate_success",
      "expected_success",
      "maximum_sample_size",
      "maximum_sample_size"
    )
  )
})

test_that("legacy calendar summaries derive each stopping reason", {
  sims <- data.frame(
    N_enrolled = c(10, 20, 30),
    stop_futility = c(0, 1, 0),
    stop_expected_success = c(1, 0, 0),
    accrual_stop_time = 1:3,
    analysis_ready_time = 4:6,
    planned_completion_time = 7:9,
    followup_person_time = 10:12,
    peak_active_followup = 1:3
  )

  expect_warning(
    out <- summarise_calendar_time(sims),
    "return_trace = TRUE"
  )
  duration <- out$trial_duration
  duration <- duration[duration$stopping_reason != "overall", ]

  expect_identical(
    duration$stopping_reason,
    c("expected_success", "futility", "maximum_sample_size")
  )
  expect_identical(duration$n_trials, c(1L, 1L, 1L))
})

test_that("calendar summaries accept factor stopping reasons", {
  sims <- calendar_test_result(include_traces = FALSE)$sims
  sims$stop_immediate_success <- c(TRUE, FALSE, FALSE)
  sims$stopping_reason <- factor(c(
    NA_character_,
    "futility",
    "maximum_sample_size"
  ))

  expect_warning(
    out <- summarise_calendar_time(sims),
    "return_trace = TRUE"
  )

  expect_equal(
    out$trial_duration$stopping_reason,
    c("immediate_success", "futility", "maximum_sample_size", "overall")
  )
})

test_that("calendar summaries retain named scenarios", {
  result <- calendar_test_result()
  out <- summarise_calendar_time(list(target = result, null = result))

  expect_setequal(unique(out$trial_duration$scenario), c("target", "null"))
  expect_setequal(unique(out$interim_timing$scenario), c("target", "null"))
})

test_that("duration remains available when interim traces were not retained", {
  result <- calendar_test_result(include_traces = FALSE)

  expect_warning(
    out <- summarise_calendar_time(result),
    "return_trace = TRUE"
  )
  expect_equal(nrow(out$trial_duration), 4L)
  expect_equal(nrow(out$interim_timing), 0L)
})

test_that("present but empty traces represent a design without interim looks", {
  result <- calendar_test_result()
  result$traces <- empty_calendar_trace()[
    0,
    setdiff(
      names(empty_calendar_trace()),
      c(".calendar_n_requested", ".calendar_n_analyzed")
    )
  ]

  expect_no_warning(out <- summarise_calendar_time(result))
  expect_equal(nrow(out$interim_timing), 0L)
})

test_that("calendar summary print is compact", {
  out <- summarise_calendar_time(calendar_test_result())

  expect_output(
    print(out),
    "Trial duration and follow-up burden"
  )
  expect_output(print(out), "median \\[P10-P90\\]")
  expect_output(print(out), "Interim-look timing")
})

test_that("old simulation results fail with a rerun instruction", {
  old <- data.frame(
    N_enrolled = 100,
    stop_futility = 0,
    stop_expected_success = 0
  )

  expect_error(
    summarise_calendar_time(old),
    "Rerun the simulations"
  )
})
