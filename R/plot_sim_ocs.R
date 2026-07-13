#' @title Plot operating characteristics across simulation scenarios
#'
#' @description Draws operating-characteristic curves across a series of true
#'   treatment-effect scenarios. The first panel shows final success and
#'   stopping probabilities. The second panel shows mean enrolled sample size.
#'
#' @param x A data frame returned by [summarise_sims()], with one row per
#'   simulation scenario.
#' @param effect Numeric treatment-effect values corresponding to the rows of
#'   `x`, or a single character string naming a numeric column in `x`.
#' @param xlab Character label for the treatment-effect axis.
#'
#' @return `x`, invisibly.
#'
#' @examples
#' operating_characteristics <- data.frame(
#'   scenario = c("null", "small", "target"),
#'   effect = c(1, 0.85, 0.7),
#'   power = c(0.025, 0.55, 0.9),
#'   stop_success = c(0.01, 0.35, 0.75),
#'   stop_futility = c(0.7, 0.25, 0.05),
#'   stop_max_N = c(0.29, 0.4, 0.2),
#'   mean_N = c(180, 250, 210),
#'   sd_N = c(45, 60, 55),
#'   stop_and_fail = c(0.001, 0.01, 0.02)
#' )
#' plot_sim_ocs(
#'   operating_characteristics,
#'   effect = "effect",
#'   xlab = "True hazard ratio"
#' )
#'
#' @export
plot_sim_ocs <- function(x, effect, xlab = "True treatment effect") {
  probability_columns <- c(
    "power", "stop_success", "stop_futility", "stop_max_N"
  )
  required <- c(probability_columns, "mean_N")
  if (!is.data.frame(x) || !all(required %in% names(x))) {
    stop(
      "'x' must be a simulation summary containing power, stopping ",
      "probabilities, and mean_N"
    )
  }
  if (nrow(x) == 0) {
    stop("'x' must contain at least one simulation scenario")
  }

  if (is.character(effect) && length(effect) == 1L) {
    if (!effect %in% names(x)) {
      stop("the effect column was not found in 'x'")
    }
    effect_values <- x[[effect]]
  } else {
    effect_values <- effect
  }
  if (
    !is.numeric(effect_values) ||
      length(effect_values) != nrow(x) ||
      any(!is.finite(effect_values))
  ) {
    stop("'effect' must supply one finite numeric value per scenario")
  }
  if (!is.character(xlab) || length(xlab) != 1L || is.na(xlab)) {
    stop("'xlab' must be a single character value")
  }

  valid_numeric <- vapply(
    x[required],
    function(value) is.numeric(value) && all(is.finite(value)),
    logical(1)
  )
  if (!all(valid_numeric)) {
    stop("operating characteristics must be finite numeric values")
  }
  probabilities <- as.matrix(x[probability_columns])
  if (any(probabilities < 0 | probabilities > 1)) {
    stop("probability operating characteristics must lie between 0 and 1")
  }
  if (any(x$mean_N < 0)) {
    stop("mean enrolled sample sizes must be non-negative")
  }

  scenario_order <- order(effect_values)
  effect_values <- effect_values[scenario_order]
  probabilities <- probabilities[scenario_order, , drop = FALSE]
  mean_N <- x$mean_N[scenario_order]

  probability_labels <- c(
    "Final success",
    "Stop for expected success",
    "Stop for futility",
    "Maximum sample size"
  )
  colours <- c("#0072B2", "#009E73", "#D55E00", "#999999")
  point_characters <- c(16, 17, 15, 18)
  line_types <- c(1, 2, 2, 3)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  graphics::matplot(
    effect_values,
    probabilities,
    type = "b",
    ylim = c(0, 1),
    xlab = xlab,
    ylab = "Probability",
    main = "Success and stopping probabilities",
    col = colours,
    pch = point_characters,
    lty = line_types
  )
  graphics::legend(
    "topright",
    legend = probability_labels,
    col = colours,
    pch = point_characters,
    lty = line_types,
    bty = "n"
  )

  y_max <- max(mean_N)
  if (y_max == 0) {
    y_max <- 1
  } else {
    y_max <- 1.08 * y_max
  }
  graphics::plot(
    effect_values,
    mean_N,
    type = "b",
    ylim = c(0, y_max),
    xlab = xlab,
    ylab = "Mean enrolled sample size",
    main = "Expected sample size",
    col = "#0072B2",
    pch = 16
  )

  invisible(x)
}
