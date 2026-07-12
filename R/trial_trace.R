#' @title Create an interim decision trace
#'
#' @description Binds trace rows produced during an adaptive trial and supplies
#'   a stable zero-row data frame when a design has no interim looks.
#'
#' @param rows List of interim trace rows.
#'
#' @return A data frame with one row per completed interim look.
#'
#' @keywords internal
#' @noRd
new_trial_trace <- function(rows) {
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) > 0) {
    return(do.call(rbind, rows))
  }

  data.frame(
    look = integer(),
    planned_N = integer(),
    calendar_time = numeric(),
    N_enrolled = integer(),
    N_treatment = integer(),
    N_control = integer(),
    events_treatment = integer(),
    events_control = integer(),
    N_pending = integer(),
    N_not_enrolled = integer(),
    ppp_stop_now = numeric(),
    success_threshold = numeric(),
    ppp_success_at_max = numeric(),
    futility_threshold = numeric(),
    decision = character(),
    warning_count = integer(),
    warning_messages = character(),
    stringsAsFactors = FALSE
  )
}

#' @title Extract an interim decision trace
#'
#' @description Accepts either a goldilocks_trial object returned by
#'   survival_adapt with return_trace enabled or the trace data frame itself.
#'
#' @param x Trial result or decision trace.
#'
#' @return An interim decision trace data frame.
#'
#' @keywords internal
#' @noRd
get_trial_trace <- function(x) {
  if (inherits(x, "goldilocks_trial")) {
    return(x$trace)
  }
  if (is.data.frame(x)) {
    return(x)
  }
  stop(
    "'x' must be a goldilocks_trial object or an interim trace data frame"
  )
}

#' @title Print an adaptive trial trace result
#'
#' @description Prints the final trial summary and reports how many interim
#'   looks were completed when survival_adapt returns a goldilocks_trial object.
#'
#' @param x A goldilocks_trial object.
#' @param ... Additional arguments passed to print.data.frame.
#'
#' @return The input object, invisibly.
#'
#' @export
print.goldilocks_trial <- function(x, ...) {
  cat("Goldilocks adaptive trial\n")
  print(x$summary, ...)
  cat("\nInterim looks completed: ", nrow(x$trace), "\n", sep = "")
  invisible(x)
}

#' @title Summarize an interim decision path
#'
#' @description Creates a compact one-row summary of the final interim look and
#'   stopping decision. Pass the goldilocks_trial object to include final
#'   analysis information, or pass its trace element to summarize the path only.
#'
#' @param x A goldilocks_trial object or an interim trace data frame.
#'
#' @return A one-row data frame.
#'
#' @export
summarise_trial_trace <- function(x) {
  trace <- get_trial_trace(x)
  summary <- if (inherits(x, "goldilocks_trial")) x$summary else NULL

  if (nrow(trace) == 0) {
    return(data.frame(
      interim_looks_completed = 0L,
      last_look = NA_integer_,
      last_decision = "no_interim_looks",
      final_N = if (is.null(summary)) NA_integer_ else summary$N_enrolled,
      final_post_prob_ha = if (is.null(summary)) NA_real_ else summary$post_prob_ha,
      stringsAsFactors = FALSE
    ))
  }

  last <- trace[nrow(trace), , drop = FALSE]
  data.frame(
    interim_looks_completed = nrow(trace),
    last_look = last$look,
    last_decision = last$decision,
    final_N = if (is.null(summary)) NA_integer_ else summary$N_enrolled,
    final_post_prob_ha = if (is.null(summary)) NA_real_ else summary$post_prob_ha,
    ppp_stop_now = last$ppp_stop_now,
    ppp_success_at_max = last$ppp_success_at_max,
    warning_count = sum(trace$warning_count),
    stringsAsFactors = FALSE
  )
}

#' @title Plot predictive probabilities and enrollment at interim looks
#'
#' @description Draws three base-R panels showing the predictive probability of
#'   success if accrual stops now, the predictive probability of success at the
#'   maximum sample size, and enrollment and observed events by treatment arm.
#'   Thresholds and early stopping decisions are marked on the probability
#'   panels.
#'
#' @param x A goldilocks_trial object or an interim trace data frame.
#'
#' @return The trace data frame, invisibly.
#'
#' @importFrom graphics plot.new title
#' @export
plot_trial_trace <- function(x) {
  trace <- get_trial_trace(x)
  if (nrow(trace) == 0) {
    plot.new()
    title(main = "No interim looks were planned")
    return(invisible(trace))
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  # Leave room for each panel's x-axis title on compact graphics devices.
  graphics::par(mfrow = c(3, 1), mar = c(5.1, 4, 2.5, 1))

  look <- trace$look
  xlim <- if (length(look) == 1) look + c(-0.5, 0.5) else range(look)
  x_ticks <- seq.int(ceiling(xlim[1]), floor(xlim[2]))
  stop_rows <- trace$decision != "continue"

  graphics::plot(
    look,
    trace$ppp_stop_now,
    type = "b",
    ylim = c(0, 1),
    xlim = xlim,
    xaxt = "n",
    xlab = "Interim look",
    ylab = "Predictive probability",
    main = "Success if accrual stops now",
    col = "#0072B2",
    pch = 16
  )
  graphics::axis(1, at = x_ticks)
  graphics::lines(
    look,
    trace$success_threshold,
    type = "b",
    col = "#D55E00",
    pch = 1,
    lty = 2
  )
  graphics::points(
    look[stop_rows],
    trace$ppp_stop_now[stop_rows],
    col = "#009E73",
    pch = 17,
    cex = 1.2
  )
  graphics::legend(
    "bottomright",
    legend = c("Predictive probability", "Success threshold", "Stopping look"),
    col = c("#0072B2", "#D55E00", "#009E73"),
    pch = c(16, 1, 17),
    lty = c(1, 2, NA),
    bty = "n"
  )

  graphics::plot(
    look,
    trace$ppp_success_at_max,
    type = "n",
    ylim = c(0, 1),
    xlim = xlim,
    xaxt = "n",
    xlab = "Interim look",
    ylab = "Predictive probability",
    main = "Success if accrual continues to maximum sample size"
  )
  graphics::axis(1, at = x_ticks)
  if (any(is.finite(trace$ppp_success_at_max))) {
    graphics::lines(
      look,
      trace$ppp_success_at_max,
      type = "b",
      col = "#0072B2",
      pch = 16
    )
    graphics::lines(
      look,
      trace$futility_threshold,
      type = "b",
      col = "#CC79A7",
      pch = 1,
      lty = 2
    )
    graphics::legend(
      "bottomright",
      legend = c("Predictive probability", "Futility threshold"),
      col = c("#0072B2", "#CC79A7"),
      pch = c(16, 1),
      lty = c(1, 2),
      bty = "n"
    )
  } else {
    graphics::text(
      x = mean(xlim),
      y = 0.5,
      labels = "Futility monitoring disabled"
    )
  }

  y_max <- max(
    trace$N_treatment,
    trace$N_control,
    trace$events_treatment,
    trace$events_control
  )
  graphics::plot(
    look,
    trace$N_treatment,
    type = "b",
    ylim = c(0, y_max),
    xlim = xlim,
    xaxt = "n",
    xlab = "Interim look",
    ylab = "Subjects",
    main = "Enrollment and observed events by arm",
    col = "#0072B2",
    pch = 16
  )
  graphics::axis(1, at = x_ticks)
  graphics::lines(look, trace$N_control, type = "b", col = "#D55E00", pch = 16)
  graphics::lines(
    look,
    trace$events_treatment,
    type = "b",
    col = "#56B4E9",
    pch = 1,
    lty = 2
  )
  graphics::lines(
    look,
    trace$events_control,
    type = "b",
    col = "#E69F00",
    pch = 1,
    lty = 2
  )
  graphics::legend(
    "topleft",
    legend = c("Treatment enrolled", "Control enrolled", "Treatment events", "Control events"),
    col = c("#0072B2", "#D55E00", "#56B4E9", "#E69F00"),
    pch = c(16, 16, 1, 1),
    lty = c(1, 1, 2, 2),
    bty = "n"
  )

  invisible(trace)
}

#' @title Plot stopping outcomes from trial simulations
#'
#' @description Draws a bar chart of expected-success, futility, and
#'   maximum-sample-size outcomes together with a histogram of enrolled sample
#'   sizes. The input can be the sims element returned by sim_trials or the
#'   complete sim_trials result.
#'
#' @param x A simulation result data frame or the list returned by sim_trials.
#'
#' @return The simulation result data frame, invisibly.
#'
#' @export
plot_sim_stopping <- function(x) {
  sims <- if (is.list(x) && !is.data.frame(x) && "sims" %in% names(x)) {
    x$sims
  } else {
    x
  }
  required <- c("stop_expected_success", "stop_futility", "N_enrolled")
  if (!is.data.frame(sims) || !all(required %in% names(sims))) {
    stop(
      "'x' must be a simulation result with stopping indicators and N_enrolled"
    )
  }

  outcome <- ifelse(
    sims$stop_expected_success,
    "Expected success",
    ifelse(sims$stop_futility, "Futility", "Maximum sample size")
  )
  outcome <- factor(
    outcome,
    levels = c("Expected success", "Futility", "Maximum sample size")
  )
  probabilities <- prop.table(table(outcome))

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  graphics::barplot(
    probabilities,
    ylim = c(0, 1),
    ylab = "Proportion of trials",
    main = "Stopping outcomes",
    col = c("#009E73", "#D55E00", "#999999")
  )
  graphics::hist(
    sims$N_enrolled,
    xlab = "Enrolled sample size",
    main = "Final enrolled sample size",
    col = "#56B4E9",
    border = "white"
  )

  invisible(sims)
}
