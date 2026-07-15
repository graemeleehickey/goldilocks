trace_args <- function() {
  list(
    hazard_treatment = -log(0.85) / 24,
    hazard_control = -log(0.7) / 24,
    cutpoints = 0,
    N_total = 80,
    lambda = 8,
    lambda_time = 0,
    interim_look = c(40, 60),
    end_of_study = 24,
    prior = c(0.1, 0.1),
    block = 2,
    rand_ratio = c(1, 1),
    prop_loss = 0.05,
    alternative = "less",
    h0 = 0,
    Fn = c(0.05, 0.05),
    Sn = c(0.95, 0.90),
    prob_ha = 0.95,
    N_impute = 2,
    N_mcmc = 2,
    method = "bayes-surv"
  )
}

run_traced_trial <- function(...) {
  set.seed(5701)
  args <- trace_args()
  overrides <- list(...)
  for (name in names(overrides)) {
    args[[name]] <- overrides[[name]]
  }
  args$return_trace <- TRUE
  do.call(survival_adapt, args)
}

test_that("survival_adapt keeps its default data-frame return value", {
  set.seed(5701)
  out <- do.call(survival_adapt, trace_args())

  expect_s3_class(out, "data.frame")
  expect_false(inherits(out, "goldilocks_trial"))
})

test_that("survival_adapt returns an auditable interim trace on request", {
  out <- run_traced_trial()
  expected_columns <- c(
    "look", "planned_N", "calendar_time", "N_enrolled", "N_treatment",
    "N_control", "events_treatment", "events_control", "N_pending",
    "N_not_enrolled", "ppp_stop_now", "success_threshold",
    "ppp_success_at_max", "futility_threshold", "decision",
    "warning_count", "warning_messages"
  )

  expect_s3_class(out, "goldilocks_trial")
  expect_s3_class(out$summary, "data.frame")
  expect_s3_class(out$trace, "data.frame")
  expect_identical(names(out$trace), expected_columns)
  expect_gte(nrow(out$trace), 1)
  expect_lte(nrow(out$trace), 2)
  expect_true(all(out$trace$decision %in% c(
    "continue", "stop_expected_success", "stop_futility"
  )))
  expect_equal(
    out$trace$N_treatment + out$trace$N_control,
    out$trace$N_enrolled
  )

  last_decision <- out$trace$decision[nrow(out$trace)]
  if (out$summary$stop_expected_success == 1) {
    expect_identical(last_decision, "stop_expected_success")
  } else if (out$summary$stop_futility == 1) {
    expect_identical(last_decision, "stop_futility")
  }
})

test_that("traced trials are reproducible with a fixed seed", {
  expect_identical(run_traced_trial(), run_traced_trial())
})

test_that("trace marks disabled futility monitoring as unavailable", {
  out <- run_traced_trial(Fn = 0)

  expect_true(all(is.na(out$trace$ppp_success_at_max)))
  expect_true(all(is.na(out$trace$futility_threshold)))
})

test_that("trace captures warnings emitted during interim analysis", {
  args <- trace_args()
  args$cutpoints <- c(0, 12)
  args$hazard_treatment <- c(-log(0.85) / 12, -log(0.85) / 12)
  args$hazard_control <- c(-log(0.7) / 12, -log(0.7) / 12)
  args$N_impute <- 1
  args$N_mcmc <- 1
  args$empty_interval <- "propagate"

  set.seed(5701)
  out <- suppressWarnings(
    do.call(survival_adapt, c(args, list(return_trace = TRUE)))
  )

  expect_true(any(out$trace$warning_count > 0))
  expect_match(paste(out$trace$warning_messages, collapse = " "), "zero subjects")
})

test_that("traces support no-interim designs", {
  out <- run_traced_trial(interim_look = NULL)
  summary <- summarise_trial_trace(out)

  expect_equal(nrow(out$trace), 0)
  expect_equal(summary$interim_looks_completed, 0)
  expect_identical(summary$last_decision, "no_interim_looks")
})

test_that("trace summaries and plots accept trial outputs", {
  out <- run_traced_trial()
  sim_data <- data.frame(
    stop_expected_success = c(TRUE, FALSE, FALSE),
    stop_futility = c(FALSE, TRUE, FALSE),
    N_enrolled = c(40, 40, 80)
  )
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_s3_class(summarise_trial_trace(out), "data.frame")
  expect_silent(plot_trial_trace(out))
  expect_silent(plot_sim_stopping(sim_data))
})

test_that("simulation stopping plot stacks outcomes by sample size", {
  sim_data <- data.frame(
    stop_expected_success = c(TRUE, FALSE, FALSE, TRUE, FALSE),
    stop_futility = c(FALSE, TRUE, FALSE, FALSE, TRUE),
    N_enrolled = c(40, 40, 80, 80, 80)
  )
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit(grDevices::dev.off(), add = TRUE)
  captured <- new.env(parent = emptyenv())
  local_mocked_bindings(
    barplot = function(height, ...) {
      captured$height <- height
      seq_len(ncol(height))
    },
    text = function(x, y, labels, ...) {
      captured$text_x <- x
      captured$text_y <- y
      captured$labels <- labels
      captured$text_args <- list(...)
    },
    mtext = function(text, ...) {
      captured$subtitle <- text
      captured$subtitle_args <- list(...)
    },
    legend = function(...) NULL,
    .package = "graphics"
  )

  expect_silent(plot_sim_stopping(sim_data))
  expect_equal(
    unname(captured$height),
    matrix(c(0.2, 0.2, 0, 0.2, 0.2, 0.2), nrow = 3)
  )
  expect_equal(unname(captured$text_y), c(0.4, 0.6))
  expect_identical(captured$labels, c("40.0%", "60.0%"))
  expect_identical(captured$text_args$cex, 0.75)
  expect_match(captured$subtitle, "^Marginal percentage")
  expect_identical(captured$subtitle_args$adj, 0)
})

test_that("simulation stopping plot supports conditional percentages", {
  sim_data <- data.frame(
    stop_expected_success = c(TRUE, FALSE, FALSE, TRUE, FALSE),
    stop_futility = c(FALSE, TRUE, FALSE, FALSE, TRUE),
    N_enrolled = c(40, 40, 80, 80, 80)
  )
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit(grDevices::dev.off(), add = TRUE)
  captured <- new.env(parent = emptyenv())
  local_mocked_bindings(
    barplot = function(height, ...) {
      captured$height <- height
      captured$args <- list(...)
      seq_len(ncol(height))
    },
    text = function(x, y, labels, ...) {
      captured$text_y <- y
      captured$labels <- labels
    },
    mtext = function(text, ...) {
      captured$subtitle <- text
    },
    legend = function(...) NULL,
    .package = "graphics"
  )

  expect_silent(plot_sim_stopping(sim_data, type = "conditional"))
  expect_equal(
    unname(captured$height),
    matrix(c(0.2, 0.2, 0, 1 / 3, 1 / 3, 1 / 3), nrow = 3)
  )
  expect_equal(unname(captured$text_y), c(0.4, 1))
  expect_identical(captured$labels, c("40.0%", "100.0%"))
  expect_match(captured$subtitle, "still active")
})

test_that("simulation stopping plot supports cumulative percentages", {
  sim_data <- data.frame(
    stop_expected_success = c(TRUE, FALSE, FALSE, TRUE, FALSE),
    stop_futility = c(FALSE, TRUE, FALSE, FALSE, TRUE),
    N_enrolled = c(40, 40, 80, 80, 80)
  )
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit(grDevices::dev.off(), add = TRUE)
  captured <- new.env(parent = emptyenv())
  local_mocked_bindings(
    barplot = function(height, ...) {
      captured$height <- height
      captured$args <- list(...)
      seq_len(ncol(height))
    },
    text = function(x, y, labels, ...) {
      captured$labels <- labels
      captured$text_args <- list(...)
    },
    mtext = function(text, ...) {
      captured$subtitle <- text
    },
    legend = function(...) NULL,
    .package = "graphics"
  )

  expect_silent(plot_sim_stopping(sim_data, type = "cumulative"))
  expect_equal(
    unname(captured$height),
    matrix(c(0.2, 0.2, 0, 0.6, 0.4, 0.4, 0.2, 0), nrow = 4)
  )
  expect_identical(
    rownames(captured$height),
    c(
      "Expected success", "Futility", "Maximum sample size",
      "Continue to next look"
    )
  )
  expect_equal(unname(colSums(captured$height)), c(1, 1))
  expect_true(all(apply(captured$height[1:3, , drop = FALSE], 1, function(x) {
    all(diff(x) >= 0)
  })))
  expect_match(captured$subtitle, "every bar sums to 100%")
  expect_identical(captured$args$ylab, "Cumulative proportion of trials")
  expect_identical(captured$text_args$cex, 0.75)
  expect_true(all(c("20.0%", "40.0%", "60.0%") %in% captured$labels))
})

test_that("simulation stopping plot validates its type", {
  sim_data <- data.frame(
    stop_expected_success = FALSE,
    stop_futility = FALSE,
    N_enrolled = 80
  )

  expect_error(plot_sim_stopping(sim_data, type = "unknown"), "arg")
})

test_that("simulation stopping legend starts beyond the plotting region", {
  sim_data <- data.frame(
    stop_expected_success = c(TRUE, FALSE, FALSE, TRUE, FALSE),
    stop_futility = c(FALSE, TRUE, FALSE, FALSE, TRUE),
    N_enrolled = c(40, 40, 80, 80, 80)
  )
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit(grDevices::dev.off(), add = TRUE)
  captured <- new.env(parent = emptyenv())
  local_mocked_bindings(
    legend = function(x, ...) {
      captured$x <- x
      captured$plot_right <- graphics::par("usr")[2]
    },
    .package = "graphics"
  )

  expect_silent(plot_sim_stopping(sim_data, type = "cumulative"))
  expect_gt(captured$x, captured$plot_right)
})

test_that("conditional and cumulative views include trace-recorded looks", {
  simulation_result <- list(
    sims = data.frame(
      stop_expected_success = c(TRUE, FALSE, FALSE, TRUE, FALSE),
      stop_futility = c(FALSE, TRUE, FALSE, FALSE, TRUE),
      N_enrolled = c(40, 40, 80, 80, 80)
    ),
    traces = data.frame(N_enrolled = c(20, 40, 80))
  )
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit(grDevices::dev.off(), add = TRUE)
  captured <- new.env(parent = emptyenv())
  captured$heights <- list()
  local_mocked_bindings(
    barplot = function(height, ...) {
      captured$heights[[length(captured$heights) + 1L]] <- height
      seq_len(ncol(height))
    },
    text = function(...) NULL,
    mtext = function(...) NULL,
    legend = function(...) NULL,
    .package = "graphics"
  )

  graphics::plot.new()
  expect_silent(plot_sim_stopping(simulation_result, type = "conditional"))
  graphics::plot.new()
  expect_silent(plot_sim_stopping(simulation_result, type = "cumulative"))

  conditional <- captured$heights[[1]]
  cumulative <- captured$heights[[2]]
  expect_identical(colnames(conditional), c("20", "40", "80"))
  expect_equal(unname(colSums(conditional)), c(0, 0.4, 1))
  expect_identical(colnames(cumulative), c("20", "40", "80"))
  expect_equal(unname(cumulative[, 1]), c(0, 0, 0, 1))
  expect_equal(unname(colSums(cumulative)), c(1, 1, 1))
})

test_that("simulation stopping plot types handle one observed sample size", {
  sim_data <- data.frame(
    stop_expected_success = c(FALSE, FALSE),
    stop_futility = c(FALSE, FALSE),
    N_enrolled = c(80, 80)
  )
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_silent(plot_sim_stopping(sim_data, type = "marginal"))
  expect_silent(plot_sim_stopping(sim_data, type = "conditional"))
  expect_silent(plot_sim_stopping(sim_data, type = "cumulative"))
})

test_that("simulation stopping flowchart follows counts through every look", {
  skip_if_not_installed("DiagrammeR")
  sim_data <- data.frame(
    stop_expected_success = c(TRUE, FALSE, TRUE, FALSE, FALSE),
    stop_futility = c(FALSE, TRUE, FALSE, TRUE, FALSE),
    N_enrolled = c(40, 40, 80, 80, 100),
    N_max = 100
  )
  simulation_result <- list(
    sims = sim_data,
    traces = data.frame(N_enrolled = c(40, 60, 80))
  )

  flowchart <- plot_sim_stopping(simulation_result, type = "flowchart")
  dot <- flowchart$x$diagram

  expect_s3_class(flowchart, "grViz")
  expect_match(dot, "All simulated trials\\nn = 5", fixed = TRUE)
  expect_match(
    dot,
    "Look 1 (N = 40)\\nContinue enrolling\\nn = 3",
    fixed = TRUE
  )
  expect_match(
    dot,
    "Look 2 (N = 60)\\nContinue enrolling\\nn = 3",
    fixed = TRUE
  )
  expect_match(
    dot,
    "Look 3 (N = 80)\\nContinue enrolling\\nn = 1",
    fixed = TRUE
  )
  expect_match(
    dot,
    "Maximum sample size (N = 100)\\nReach maximum sample size\\nn = 1",
    fixed = TRUE
  )
  expect_match(dot, "Stop for futility\\nn = 1", fixed = TRUE)
  expect_match(dot, "Stop for early success\\nn = 1", fixed = TRUE)
})
