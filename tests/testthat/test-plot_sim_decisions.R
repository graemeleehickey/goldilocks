sim_decision_traces <- function() {
  data.frame(
    trial = rep(1:4, 2),
    look = rep(c(1, 2), each = 4),
    planned_N = rep(c(40, 60), each = 4),
    ppp_stop_now = c(0.96, 0.4, 0.2, 0.4, 0.92, 0.5, 0.1, 0.5),
    success_threshold = rep(c(0.95, 0.9), each = 4),
    ppp_success_at_max = c(0.8, 0.5, 0.02, 0.5, 0.85, 0.4, 0.01, 0.4),
    futility_threshold = 0.05,
    decision = c(
      "stop_expected_success", "continue", "stop_futility", "continue",
      "stop_expected_success", "continue", "stop_futility", "continue"
    ),
    stringsAsFactors = FALSE
  )
}

test_that("simulation decision map accepts retained traces", {
  traces <- sim_decision_traces()
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_identical(plot_sim_decisions(traces), traces)
  expect_identical(plot_sim_decisions(list(traces = traces)), traces)
})

test_that("simulation decision map handles disabled futility monitoring", {
  traces <- sim_decision_traces()
  traces$ppp_success_at_max <- NA_real_
  traces$futility_threshold <- NA_real_
  traces$decision <- "continue"
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_silent(plot_sim_decisions(traces))
})

test_that("simulation decision map handles designs without interim looks", {
  traces <- data.frame(
    look = numeric(),
    ppp_stop_now = numeric(),
    success_threshold = numeric(),
    ppp_success_at_max = numeric(),
    futility_threshold = numeric(),
    decision = character()
  )
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_silent(plot_sim_decisions(traces))
})

test_that("simulation decision map validates trace structure", {
  traces <- sim_decision_traces()
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_error(
    plot_sim_decisions(traces[, -4]),
    "simulation traces"
  )
  traces$success_threshold[1] <- 0.9
  expect_error(
    plot_sim_decisions(traces),
    "constant within each interim look"
  )
})
