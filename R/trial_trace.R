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
#' @description Draws a stacked bar chart of stopping outcomes by enrolled
#'   sample size, with colours distinguishing expected-success, futility, and
#'   maximum-sample-size outcomes. The `type` argument controls whether the
#'   function draws marginal, conditional, or cumulative bars, or a flowchart
#'   through successive interim looks. Bar-chart subtitles state the
#'   denominator used by the selected view. The input can be the `sims` element
#'   returned by [sim_trials()] or the complete `sim_trials()` result.
#'
#' @param x A simulation result data frame or the list returned by sim_trials.
#' @param type Character string specifying the percentages to plot. `"marginal"`
#'   shows the percentage of all simulated trials ending at each sample size;
#'   its bars sum to 100 percent across sample sizes. `"conditional"` shows the
#'   percentage stopping at each look among trials still active at the start of
#'   that look. `"cumulative"` shows the status of all simulated trials after
#'   each look; every bar sums to 100 percent and includes trials continuing to
#'   the next look. `"flowchart"` starts with all simulated trials and branches
#'   at each look into futility, continued enrollment, and expected-success
#'   nodes labelled with trial counts.
#'
#' @details The marginal view uses terminal sample sizes observed in
#'   `N_enrolled`. When the complete result from
#'   `sim_trials(return_trace = TRUE)` is supplied, the conditional, cumulative,
#'   and flowchart views also include sample sizes recorded in `traces`, so
#'   reached looks at which no trial stopped still appear. The flowchart
#'   requires the `N_max` column and is rendered with [DiagrammeR::grViz()].
#'
#' @return For bar-chart types, the simulation result data frame, invisibly.
#'   For `type = "flowchart"`, a `DiagrammeR` `grViz` htmlwidget.
#'
#' @export
plot_sim_stopping <- function(
  x,
  type = c("marginal", "conditional", "cumulative", "flowchart")
) {
  type <- match.arg(type)
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
  trace_sizes <- sim_stopping_trace_sizes(x)

  if (type == "flowchart") {
    return(plot_sim_stopping_flowchart(sims, trace_sizes = trace_sizes))
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
  sample_sizes <- if (type == "marginal") {
    sort(unique(sims$N_enrolled))
  } else {
    sort(unique(c(sims$N_enrolled, trace_sizes)))
  }
  sample_size <- factor(sims$N_enrolled, levels = sample_sizes)
  counts <- unclass(table(outcome, sample_size))

  if (type == "marginal") {
    probabilities <- counts / nrow(sims)
    bar_percentages <- colSums(probabilities)
    subtitle <- paste(
      "Marginal percentage of all simulated trials;",
      "bars sum to 100% across sample sizes"
    )
    ylab <- "Proportion of all simulated trials"
    xlab <- "Final enrolled sample size"
    main <- "Stopping outcomes by sample size"
  } else if (type == "conditional") {
    risk_set <- vapply(
      sample_sizes,
      function(size) sum(sims$N_enrolled >= size),
      numeric(1)
    )
    probabilities <- sweep(counts, 2, risk_set, "/")
    bar_percentages <- colSums(probabilities)
    subtitle <- paste(
      "Conditional percentage among trials still active",
      "at the start of each look"
    )
    ylab <- "Proportion of trials entering the look"
    xlab <- "Enrolled sample size at look"
    main <- "Stopping outcomes by sample size"
  } else {
    cumulative_outcomes <- counts
    for (outcome_row in seq_len(nrow(cumulative_outcomes))) {
      cumulative_outcomes[outcome_row, ] <- cumsum(counts[outcome_row, ])
    }
    cumulative_outcomes <- cumulative_outcomes / nrow(sims)
    continue <- pmax(0, 1 - colSums(cumulative_outcomes))
    probabilities <- rbind(
      cumulative_outcomes,
      "Continue to next look" = continue
    )
    subtitle <- paste(
      "Cumulative status after each look;",
      "every bar sums to 100% of simulated trials"
    )
    ylab <- "Cumulative proportion of trials"
    xlab <- "Enrolled sample size at look"
    main <- "Cumulative stopping outcomes by sample size"
  }

  bar_colours <- c("#009E73", "#D55E00", "#999999", "#56B4E9")
  percentage_cex <- 0.75
  legend_order <- rev(seq_len(nrow(probabilities)))
  legend_cex <- 0.9
  legend_width <- max(graphics::strwidth(
    rownames(probabilities)[legend_order],
    units = "inches",
    cex = legend_cex
  ))
  legend_margin <- ceiling(
    (legend_width + 0.5) / graphics::par("csi")
  )

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(
    mfrow = c(1, 1),
    mar = c(5, 4, 4, max(6, legend_margin))
  )
  bar_midpoints <- graphics::barplot(
    probabilities,
    ylim = c(0, 1.08),
    xlab = xlab,
    ylab = ylab,
    main = main,
    col = bar_colours,
    border = NA
  )
  graphics::mtext(
    subtitle,
    side = 3,
    line = 0.2,
    adj = 0,
    cex = 0.85
  )
  plot_region <- graphics::par("usr")
  graphics::legend(
    x = plot_region[2] + 0.04 * diff(plot_region[1:2]),
    y = plot_region[4],
    legend = rownames(probabilities)[legend_order],
    fill = bar_colours[legend_order],
    border = NA,
    bty = "n",
    xjust = 0,
    yjust = 1,
    xpd = NA,
    cex = legend_cex
  )
  if (type == "cumulative") {
    segment_midpoints <- apply(probabilities, 2, cumsum) - probabilities / 2
    segment_labels <- ifelse(
      probabilities > 0,
      sprintf("%.1f%%", 100 * probabilities),
      ""
    )
    graphics::text(
      rep(bar_midpoints, each = nrow(probabilities)),
      as.vector(segment_midpoints),
      labels = as.vector(segment_labels),
      cex = percentage_cex
    )
  } else {
    graphics::text(
      bar_midpoints,
      bar_percentages,
      labels = sprintf("%.1f%%", 100 * bar_percentages),
      pos = 3,
      offset = 0.25,
      cex = percentage_cex
    )
  }

  invisible(sims)
}

#' Extract sample sizes from retained simulation traces
#'
#' @param simulation_result Original input supplied to `plot_sim_stopping()`.
#'
#' @return A numeric vector of trace-recorded enrolled sample sizes.
#'
#' @noRd
sim_stopping_trace_sizes <- function(simulation_result) {
  if (
    is.list(simulation_result) &&
      !is.data.frame(simulation_result) &&
      "traces" %in% names(simulation_result) &&
      is.data.frame(simulation_result$traces) &&
      "N_enrolled" %in% names(simulation_result$traces)
  ) {
    sizes <- simulation_result$traces$N_enrolled
    return(sizes[is.finite(sizes)])
  }

  numeric()
}

#' Draw a stopping flowchart for simulated trials
#'
#' @param sims Simulation summary data frame.
#' @param trace_sizes Sample sizes from retained simulation traces.
#'
#' @return A `DiagrammeR` `grViz` htmlwidget.
#'
#' @noRd
plot_sim_stopping_flowchart <- function(sims, trace_sizes) {
  if (!requireNamespace("DiagrammeR", quietly = TRUE)) {
    stop(
      "Package 'DiagrammeR' is required for type = 'flowchart'. ",
      "Install it with install.packages('DiagrammeR')."
    )
  }

  if (!"N_max" %in% names(sims)) {
    stop("type = 'flowchart' requires the N_max simulation column")
  }
  planned_max <- unique(sims$N_max)
  if (length(planned_max) != 1L || !is.finite(planned_max)) {
    stop("flowcharts require simulations from a single maximum sample size")
  }

  sample_sizes <- sort(unique(c(
    sims$N_enrolled,
    trace_sizes,
    planned_max
  )))
  sample_sizes <- sample_sizes[is.finite(sample_sizes)]
  if (length(sample_sizes) == 0L) {
    stop("flowcharts require at least one observed sample size")
  }
  if (any(sample_sizes > planned_max)) {
    stop("interim sample sizes cannot exceed N_max")
  }

  reached_max <- !sims$stop_expected_success & !sims$stop_futility
  if (any(sims$N_enrolled[reached_max] != planned_max)) {
    stop("maximum-sample-size outcomes must share the same N_max")
  }
  stopped_early <- sims$stop_expected_success | sims$stop_futility
  if (any(sims$N_enrolled[stopped_early] >= planned_max)) {
    stop("early-stopping outcomes must occur before N_max")
  }
  interim_sizes <- sample_sizes[sample_sizes < planned_max]

  node_statements <- sprintf(
    paste0(
      "total [label=\"All simulated trials\\nn = %d\", ",
      "fillcolor=\"#F2F2F2\", color=\"#666666\", penwidth=1.4]"
    ),
    nrow(sims)
  )
  edge_statements <- character()
  rank_statements <- character()
  previous_continue <- "total"

  for (look_index in seq_along(interim_sizes)) {
    sample_size <- interim_sizes[look_index]
    size_label <- format(sample_size, trim = TRUE, scientific = FALSE)
    look_label <- paste("Look", look_index)
    stopped_futility <- sum(
      sims$N_enrolled == sample_size & sims$stop_futility
    )
    stopped_success <- sum(
      sims$N_enrolled == sample_size & sims$stop_expected_success
    )
    center_count <- sum(sims$N_enrolled > sample_size)

    futility_id <- paste0("futility_", look_index)
    continue_id <- paste0("continue_", look_index)
    success_id <- paste0("success_", look_index)
    node_statements <- c(
      node_statements,
      sprintf(
        paste0(
          "%s [label=\"%s (N = %s)\\nStop for futility\\nn = %d\", ",
          "fillcolor=\"#F6DED4\", color=\"#D55E00\"]"
        ),
        futility_id,
        look_label,
        size_label,
        stopped_futility
      ),
      sprintf(
        paste0(
          "%s [label=\"%s (N = %s)\\n%s\\nn = %d\", ",
          "fillcolor=\"%s\", color=\"%s\"]"
        ),
        continue_id,
        look_label,
        size_label,
        "Continue enrolling",
        center_count,
        "#D8EEF8",
        "#56B4E9"
      ),
      sprintf(
        paste0(
          "%s [label=\"%s (N = %s)\\nStop for early success\\nn = %d\", ",
          "fillcolor=\"#D5F1E7\", color=\"#009E73\"]"
        ),
        success_id,
        look_label,
        size_label,
        stopped_success
      )
    )
    edge_statements <- c(
      edge_statements,
      sprintf(
        "%s -> {%s %s %s}",
        previous_continue,
        futility_id,
        continue_id,
        success_id
      ),
      sprintf("%s -> %s [style=invis]", futility_id, continue_id),
      sprintf("%s -> %s [style=invis]", continue_id, success_id)
    )
    rank_statements <- c(
      rank_statements,
      sprintf(
        "{rank=same; %s; %s; %s}",
        futility_id,
        continue_id,
        success_id
      )
    )
    previous_continue <- continue_id
  }

  max_label <- format(planned_max, trim = TRUE, scientific = FALSE)
  node_statements <- c(
    node_statements,
    sprintf(
      paste0(
        "maximum [label=\"Maximum sample size (N = %s)\\n",
        "Reach maximum sample size\\nn = %d\", ",
        "fillcolor=\"#E6E6E6\", color=\"#777777\", penwidth=1.4]"
      ),
      max_label,
      sum(reached_max)
    )
  )
  edge_statements <- c(
    edge_statements,
    sprintf("%s -> maximum", previous_continue)
  )

  dot <- paste(
    "digraph stopping_flow {",
    paste0(
      "graph [rankdir=TB, bgcolor=\"transparent\", ranksep=0.65, ",
      "nodesep=0.35, pad=0.2]"
    ),
    paste0(
      "node [shape=box, style=\"rounded,filled\", fontname=\"Helvetica\", ",
      "fontsize=11, fontcolor=\"#222222\", margin=\"0.15,0.10\"]"
    ),
    paste0(
      "edge [color=\"#777777\", penwidth=1.2, arrowsize=0.7, ",
      "fontname=\"Helvetica\"]"
    ),
    paste(node_statements, collapse = "\n"),
    paste(edge_statements, collapse = "\n"),
    paste(rank_statements, collapse = "\n"),
    "}",
    sep = "\n"
  )

  DiagrammeR::grViz(dot)
}
