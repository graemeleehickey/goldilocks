#' @title Plot predictive-probability decision maps
#'
#' @description Draws one decision map per interim look across simulated
#'   trials. The horizontal axis is the predictive probability of success if
#'   enrollment continues to the maximum sample size, and the vertical axis is
#'   the predictive probability of success if enrollment stops at the current
#'   sample size. Shaded regions and dashed lines show the futility,
#'   continuation, stopping-accrual-for-expected-success, and immediate-success
#'   rules. Identical points are aggregated; point size indicates their
#'   frequency.
#'
#' @param x A required result returned by [sim_trials()] with
#'   `return_trace = TRUE`, or the corresponding `traces` data frame. The trace
#'   must contain interim-look numbers, predictive probabilities, decision
#'   thresholds, and decisions.
#'
#' @details The upper shaded region represents
#'   `ppp_stop_now > immediate_success_threshold`; the middle region represents
#'   `success_threshold < ppp_stop_now <= immediate_success_threshold`.
#'   Futility applies below or at the expected-success threshold when the
#'   maximum-sample-size prediction is strictly below its threshold. A legacy
#'   trace without `immediate_success_threshold` is treated as if that threshold
#'   were one, which disables immediate success.
#'
#' @return The simulation traces, invisibly.
#'
#' @examples
#' traces <- data.frame(
#'   trial = 1:6,
#'   look = rep(c(1, 2), each = 3),
#'   planned_N = rep(c(40, 60), each = 3),
#'   ppp_stop_now = c(0.96, 0.4, 0.2, 0.92, 0.5, 0.1),
#'   success_threshold = rep(c(0.95, 0.9), each = 3),
#'   immediate_success_threshold = 0.99,
#'   ppp_success_at_max = c(0.8, 0.5, 0.02, 0.85, 0.4, 0.01),
#'   futility_threshold = 0.05,
#'   decision = c(
#'     "stop_expected_success", "continue", "stop_futility",
#'     "stop_expected_success", "continue", "stop_futility"
#'   )
#' )
#' plot_sim_decisions(traces)
#'
#' @export
plot_sim_decisions <- function(x) {
  traces <- if (is.list(x) && !is.data.frame(x) && "traces" %in% names(x)) {
    x$traces
  } else {
    x
  }
  required <- c(
    "look",
    "ppp_stop_now",
    "success_threshold",
    "ppp_success_at_max",
    "futility_threshold",
    "decision"
  )
  if (!is.data.frame(traces) || !all(required %in% names(traces))) {
    stop(
      "'x' must contain simulation traces with predictive probabilities, ",
      "thresholds, and decisions"
    )
  }
  if (nrow(traces) == 0) {
    graphics::plot.new()
    graphics::title(main = "No interim decisions were simulated")
    return(invisible(traces))
  }
  returned_traces <- traces
  if (!"immediate_success_threshold" %in% names(traces)) {
    traces$immediate_success_threshold <- rep.int(1, nrow(traces))
  }

  finite_probability <- function(value) {
    is.numeric(value) &&
      all(is.finite(value)) &&
      all(value >= 0 & value <= 1)
  }
  optional_probability <- function(value) {
    is.numeric(value) &&
      all(is.na(value) | is.finite(value)) &&
      all(value[!is.na(value)] >= 0 & value[!is.na(value)] <= 1)
  }
  if (!is.numeric(traces$look) || any(!is.finite(traces$look))) {
    stop("trace look numbers must be finite numeric values")
  }
  if (
    !finite_probability(traces$ppp_stop_now) ||
      !finite_probability(traces$success_threshold) ||
      !finite_probability(traces$immediate_success_threshold) ||
      !optional_probability(traces$ppp_success_at_max) ||
      !optional_probability(traces$futility_threshold)
  ) {
    stop("predictive probabilities and thresholds must lie between 0 and 1")
  }
  allowed_decisions <- c(
    "continue",
    "stop_immediate_success",
    "stop_expected_success",
    "stop_futility"
  )
  if (
    !is.character(traces$decision) ||
      any(is.na(traces$decision)) ||
      !all(traces$decision %in% allowed_decisions)
  ) {
    stop("trace decisions contain an unsupported value")
  }
  futility_available <- is.finite(traces$ppp_success_at_max) &
    is.finite(traces$futility_threshold)
  if (
    any(xor(
      is.finite(traces$ppp_success_at_max),
      is.finite(traces$futility_threshold)
    ))
  ) {
    stop("maximum-sample-size predictions and futility thresholds must align")
  }

  looks <- sort(unique(traces$look))
  panel_columns <- ceiling(sqrt(length(looks)))
  panel_rows <- ceiling(length(looks) / panel_columns)
  decision_labels <- c(
    "Continue",
    "Immediate success",
    "Stop accrual for expected success",
    "Futility"
  )
  decision_colours <- c("#0072B2", "#CC79A7", "#009E73", "#D55E00")
  names(decision_colours) <- allowed_decisions

  for (look in looks) {
    rows <- traces[traces$look == look, , drop = FALSE]
    if (length(unique(rows$success_threshold)) != 1L) {
      stop("success thresholds must be constant within each interim look")
    }
    if (length(unique(rows$immediate_success_threshold)) != 1L) {
      stop(
        "immediate-success thresholds must be constant within each interim look"
      )
    }
    if (any(rows$success_threshold > rows$immediate_success_threshold)) {
      stop("immediate-success thresholds cannot be below success thresholds")
    }
    look_has_futility <- futility_available[traces$look == look]
    if (any(look_has_futility) && !all(look_has_futility)) {
      stop("futility monitoring must be consistent within each interim look")
    }
    if (
      any(look_has_futility) &&
        length(unique(rows$futility_threshold)) != 1L
    ) {
      stop("futility thresholds must be constant within each interim look")
    }
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(
    mfrow = c(panel_rows, panel_columns),
    mar = c(4, 4, 3, 1)
  )

  for (i in seq_along(looks)) {
    look <- looks[i]
    rows <- traces[traces$look == look, , drop = FALSE]
    success_threshold <- unique(rows$success_threshold)
    immediate_success_threshold <- unique(
      rows$immediate_success_threshold
    )
    look_has_futility <- futility_available[traces$look == look]

    planned_N <- if (
      "planned_N" %in% names(rows) && is.numeric(rows$planned_N)
    ) {
      unique(rows$planned_N[is.finite(rows$planned_N)])
    } else {
      numeric()
    }
    panel_title <- if (length(planned_N) == 1L) {
      sprintf("Interim look %s (N = %s)", look, planned_N)
    } else {
      sprintf("Interim look %s", look)
    }

    futility_enabled <- any(look_has_futility)
    futility_threshold <- if (futility_enabled) {
      unique(rows$futility_threshold)
    } else {
      NA_real_
    }

    graphics::plot(
      NA_real_,
      NA_real_,
      type = "n",
      xlim = c(0, 1),
      ylim = c(0, 1),
      xlab = if (futility_enabled) {
        "Success probability at maximum sample size"
      } else {
        "Futility disabled (horizontal position arbitrary)"
      },
      ylab = "Success probability if enrollment stops now",
      main = panel_title
    )
    graphics::rect(
      0,
      0,
      1,
      1,
      col = grDevices::adjustcolor("#0072B2", alpha.f = 0.06),
      border = NA
    )
    if (futility_enabled) {
      graphics::rect(
        0,
        0,
        futility_threshold,
        success_threshold,
        col = grDevices::adjustcolor("#D55E00", alpha.f = 0.14),
        border = NA
      )
    }
    graphics::rect(
      0,
      success_threshold,
      1,
      immediate_success_threshold,
      col = grDevices::adjustcolor("#009E73", alpha.f = 0.14),
      border = NA
    )
    graphics::rect(
      0,
      immediate_success_threshold,
      1,
      1,
      col = grDevices::adjustcolor("#CC79A7", alpha.f = 0.14),
      border = NA
    )
    graphics::abline(
      h = success_threshold,
      col = "#009E73",
      lty = 2
    )
    graphics::abline(
      h = immediate_success_threshold,
      col = "#CC79A7",
      lty = 3
    )
    if (futility_enabled) {
      graphics::abline(
        v = futility_threshold,
        col = "#D55E00",
        lty = 2
      )
    }

    point_x <- if (futility_enabled) {
      rows$ppp_success_at_max
    } else {
      rep.int(0.5, nrow(rows))
    }
    point_data <- stats::aggregate(
      list(count = rep.int(1L, nrow(rows))),
      by = list(
        ppp_success_at_max = point_x,
        ppp_stop_now = rows$ppp_stop_now,
        decision = rows$decision
      ),
      FUN = sum
    )
    point_size <- 0.65 + 1.75 * sqrt(point_data$count / max(point_data$count))
    graphics::points(
      point_data$ppp_success_at_max,
      point_data$ppp_stop_now,
      pch = 16,
      cex = point_size,
      col = grDevices::adjustcolor(
        decision_colours[point_data$decision],
        alpha.f = 0.7
      )
    )
    if (i == 1L) {
      graphics::legend(
        "topright",
        legend = decision_labels,
        col = decision_colours,
        pch = 16,
        bty = "n",
        cex = 0.8
      )
    }
  }

  invisible(returned_traces)
}
