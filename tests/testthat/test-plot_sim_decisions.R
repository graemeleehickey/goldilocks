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
      "stop_expected_success",
      "continue",
      "stop_futility",
      "continue",
      "stop_expected_success",
      "continue",
      "stop_futility",
      "continue"
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

test_that("simulation decision map renders four exclusive regions", {
  traces <- data.frame(
    trial = 1:4,
    look = 1,
    planned_N = 40,
    ppp_stop_now = c(0.99, 0.95, 0.90, 0.20),
    success_threshold = 0.90,
    immediate_success_threshold = 0.98,
    ppp_success_at_max = c(0.80, 0.80, 0.01, 0.50),
    futility_threshold = 0.05,
    decision = c(
      "stop_immediate_success",
      "stop_expected_success",
      "stop_futility",
      "continue"
    )
  )
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit(grDevices::dev.off(), add = TRUE)
  captured <- new.env(parent = emptyenv())
  captured$rectangles <- list()
  captured$boundaries <- list()
  local_mocked_bindings(
    rect = function(xleft, ybottom, xright, ytop, ...) {
      captured$rectangles[[length(captured$rectangles) + 1L]] <-
        c(xleft, ybottom, xright, ytop)
    },
    abline = function(...) {
      captured$boundaries[[length(captured$boundaries) + 1L]] <- list(...)
    },
    .package = "graphics"
  )

  expect_identical(plot_sim_decisions(traces), traces)
  expect_equal(captured$rectangles[[2]], c(0, 0, 0.05, 0.90))
  expect_equal(captured$rectangles[[3]], c(0, 0.90, 1, 0.98))
  expect_equal(captured$rectangles[[4]], c(0, 0.98, 1, 1))
  expect_equal(captured$boundaries[[1]]$h, 0.90)
  expect_equal(captured$boundaries[[2]]$h, 0.98)
  expect_equal(captured$boundaries[[3]]$v, 0.05)
})

test_that("simulation decision map tolerates a disabled upper boundary", {
  traces <- sim_decision_traces()
  traces$immediate_success_threshold <- 1
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_silent(plot_sim_decisions(traces))
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

  traces <- sim_decision_traces()
  traces$immediate_success_threshold <- traces$success_threshold - 0.01
  expect_error(
    plot_sim_decisions(traces),
    "cannot be below success thresholds"
  )
})
