#' @title Summarize a finite Monte Carlo probability estimate
#'
#' @description Computes the point estimate, Monte Carlo standard error, and
#'   exact one-sided Clopper-Pearson bounds for a binary Monte Carlo estimand.
#'   Threshold crossing uses the point estimate; the bounds are diagnostic and
#'   do not alter the decision.
#'
#' @noRd
monte_carlo_probability_summary <- function(
  successes,
  draws,
  threshold,
  direction = c("greater", "less"),
  confidence = 0.95
) {
  direction <- match.arg(direction)
  validate_nonnegative_integer_scalar(successes, "successes")
  validate_positive_integer_scalar(draws, "draws")
  if (successes > draws) {
    stop("'successes' must not exceed 'draws'")
  }
  validate_single_probability(threshold, "threshold")
  validate_single_probability(confidence, "confidence", upper_open = TRUE)
  if (confidence <= 0.5) {
    stop("'confidence' must be greater than 0.5 and less than 1")
  }

  estimate <- successes / draws
  mcse <- sqrt(estimate * (1 - estimate) / draws)
  alpha <- 1 - confidence
  lower <- if (successes == 0L) {
    0
  } else {
    stats::qbeta(alpha, successes, draws - successes + 1L)
  }
  upper <- if (successes == draws) {
    1
  } else {
    stats::qbeta(confidence, successes + 1L, draws - successes)
  }

  point_crossed <- if (direction == "greater") {
    estimate > threshold
  } else {
    estimate < threshold
  }
  bound_crossed <- if (direction == "greater") {
    lower > threshold
  } else {
    upper < threshold
  }
  reason <- if (point_crossed) {
    paste0("estimate_", direction, "_threshold")
  } else {
    paste0("estimate_not_", direction, "_threshold")
  }

  list(
    estimate = estimate,
    mcse = mcse,
    lower = lower,
    upper = upper,
    successes = as.integer(successes),
    draws = as.integer(draws),
    threshold = threshold,
    direction = direction,
    confidence = confidence,
    point_crossed = point_crossed,
    bound_crossed = bound_crossed,
    crossed = point_crossed,
    reason = reason
  )
}

#' @title Attach binary Monte Carlo counts to an analysis result
#'
#' @noRd
set_analysis_mc_counts <- function(result, successes, draws) {
  attr(result, "mc_counts") <- list(
    successes = as.integer(successes),
    draws = as.integer(draws)
  )
  result
}

#' @title Classify one completed-data analysis
#'
#' @description All analyses use their point result. For analyses based on
#'   posterior Monte Carlo draws, exact bounds are retained as diagnostics but
#'   do not alter the completed-dataset classification.
#'
#' @noRd
classify_completed_analysis <- function(
  analysis,
  prob_ha,
  mc_conf_level
) {
  mc_counts <- attr(analysis, "mc_counts", exact = TRUE)
  if (is.null(mc_counts)) {
    crossed <- analysis$success > prob_ha
    return(list(
      crossed = crossed,
      uncertain = FALSE,
      estimate = analysis$success,
      lower = NA_real_,
      upper = NA_real_,
      draws = NA_integer_
    ))
  }

  summary <- monte_carlo_probability_summary(
    successes = mc_counts$successes,
    draws = mc_counts$draws,
    threshold = prob_ha,
    direction = "greater",
    confidence = mc_conf_level
  )
  list(
    crossed = summary$point_crossed,
    uncertain = summary$point_crossed && !summary$bound_crossed,
    estimate = summary$estimate,
    lower = summary$lower,
    upper = summary$upper,
    draws = summary$draws
  )
}
