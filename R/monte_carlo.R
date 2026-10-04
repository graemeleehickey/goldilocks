#' @title Summarize a finite Monte Carlo probability estimate
#'
#' @description Computes the point estimate, Monte Carlo standard error, and
#'   exact one-sided Clopper-Pearson bounds for a binary Monte Carlo estimand.
#'   Threshold crossing uses the point estimate; the bounds are diagnostic and
#'   do not alter the decision.
#'
#' @param successes A non-negative integer giving the number of successful Monte
#'   Carlo draws.
#' @param draws A positive integer giving the total number of Monte Carlo draws.
#' @param threshold A single numeric probability in `[0, 1]` against which the
#'   estimated probability is compared.
#' @param direction A character string specifying whether the estimate must be
#'   `"greater"` than or `"less"` than `threshold`. The default is `"greater"`.
#' @param confidence A single numeric probability strictly between `0.5` and
#'   `1`, giving the confidence level for the diagnostic bound. The default is
#'   `0.95`.
#'
#' @return A list containing the probability estimate, Monte Carlo standard
#'   error, one-sided bounds, draw counts, and threshold-crossing indicators.
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
#' @param result A completed-data analysis result.
#' @param successes A non-negative integer giving the number of posterior draws
#'   that satisfy the alternative hypothesis.
#' @param draws A positive integer giving the total number of posterior draws.
#'
#' @return `result` with the counts stored in its `mc_counts` attribute.
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
#' @param analysis A list returned by `analyse_data()`.
#' @param prob_ha A single numeric probability in `[0, 1]` defining success.
#' @param mc_conf_level A single numeric probability strictly between `0.5` and
#'   `1`, giving the confidence level for diagnostic bounds.
#'
#' @return A list containing the point-estimate classification and any
#'   diagnostic Monte Carlo bounds.
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

# Retain only the success scores and posterior counts needed for recalibration.
pack_analysis_scores <- function(analyses) {
  errors <- vapply(analyses, inherits, logical(1), what = "error")
  if (any(errors)) {
    error <- analyses[[which(errors)[1L]]]
    return(structure(
      numeric(),
      failure = list(
        error_class = class(error)[1L],
        message = conditionMessage(error)
      )
    ))
  }
  scores <- vapply(analyses, `[[`, numeric(1), "success")
  if (
    !length(scores) || any(!is.finite(scores)) || any(scores < 0 | scores > 1)
  ) {
    return(structure(
      numeric(),
      failure = list(
        error_class = "goldilocks_invalid_predictive_score",
        message = "Predictive analysis did not return finite probabilities in [0, 1]"
      )
    ))
  }
  counts <- lapply(analyses, attr, which = "mc_counts", exact = TRUE)
  has_counts <- !vapply(counts, is.null, logical(1))
  if (any(has_counts)) {
    if (!all(has_counts)) {
      stop("Inconsistent posterior counts in predictive analyses")
    }
    draws <- vapply(counts, function(x) as.integer(x$draws), integer(1))
    attr(scores, "mc_counts") <- list(
      successes = vapply(
        counts,
        function(x) as.integer(x$successes),
        integer(1)
      ),
      draws = if (length(unique(draws)) == 1L) draws[1L] else draws
    )
  }
  scores
}

classify_analysis_scores <- function(scores, prob_ha, mc_conf_level) {
  if (
    !length(scores) || any(!is.finite(scores)) || any(scores < 0 | scores > 1)
  ) {
    stop("Predictive scores must be complete, finite probabilities")
  }
  crossed <- scores > prob_ha
  counts <- attr(scores, "mc_counts", exact = TRUE)
  uncertain <- if (is.null(counts)) {
    0L
  } else {
    lower <- numeric(length(scores))
    positive <- counts$successes > 0L
    draws <- rep_len(counts$draws, length(scores))
    lower[positive] <- stats::qbeta(
      1 - mc_conf_level,
      counts$successes[positive],
      draws[positive] - counts$successes[positive] + 1L
    )
    sum(crossed & lower <= prob_ha)
  }
  list(successes = sum(crossed), uncertain = uncertain)
}

# Banks isolate current- and maximum-cohort failures; ordinary analyses still stop.
capture_predictive_analysis <- function(isolate_errors, expr) {
  if (isolate_errors) tryCatch(expr, error = identity) else expr
}
