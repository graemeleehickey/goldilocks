#' Convert a normal-test statistic to the package success scale
#'
#' @param statistic A numeric normal-score statistic.
#' @param alternative The direction of the alternative hypothesis.
#'
#' @return One minus the normal-approximation P-value.
#'
#' @keywords internal
#' @noRd
normal_test_success <- function(statistic, alternative) {
  switch(
    alternative,
    "less" = 1 - pnorm(statistic),
    "greater" = pnorm(statistic),
    "two.sided" = 1 - 2 * pnorm(-abs(statistic))
  )
}

#' @title Pool scalar treatment effects using Rubin's rules
#'
#' @description Combines scalar treatment-effect estimates and their
#'   within-imputation variances, then evaluates the pooled estimate against
#'   its null value using Rubin's large-sample degrees of freedom.
#'
#' @param estimates A numeric vector of treatment-effect estimates, one per
#'   imputed data set. At least two values are required.
#' @param variances A numeric vector of corresponding finite, non-negative
#'   within-imputation variances.
#' @inheritParams survival_adapt
#'
#' @return A list containing the pooled success score, effect estimate, standard
#'   error, and degrees of freedom.
#'
#' @importFrom stats pt var
#' @noRd
pool_rubin_scalar <- function(estimates, variances, alternative, h0) {
  m <- length(estimates)
  if (m < 2 || length(variances) != m) {
    stop("Rubin pooling requires at least two paired estimates and variances")
  }
  if (
    anyNA(estimates) ||
      anyNA(variances) ||
      any(!is.finite(estimates)) ||
      any(!is.finite(variances)) ||
      any(variances < 0)
  ) {
    stop(
      "Rubin pooling requires finite estimates and non-negative variances"
    )
  }

  estimate <- mean(estimates)
  within_variance <- mean(variances)
  between_variance <- var(estimates)
  total_variance <- within_variance + (1 + 1 / m) * between_variance
  if (!is.finite(total_variance) || total_variance < 0) {
    stop("Rubin pooling requires a finite non-negative total variance")
  }
  if (total_variance == 0) {
    difference_from_null <- estimate - h0
    statistic <- if (difference_from_null == 0) {
      0
    } else {
      sign(difference_from_null) * Inf
    }
    return(list(
      success = normal_test_success(statistic, alternative),
      estimate = estimate,
      std_error = 0,
      degrees_freedom = Inf
    ))
  }

  relative_increase <- if (within_variance == 0) {
    Inf
  } else {
    (1 + 1 / m) * between_variance / within_variance
  }
  degrees_freedom <- if (relative_increase == 0) {
    Inf
  } else {
    (m - 1) * (1 + 1 / relative_increase)^2
  }

  statistic <- (estimate - h0) / sqrt(total_variance)
  success <- switch(
    alternative,
    "less" = 1 - pt(statistic, df = degrees_freedom),
    "greater" = pt(statistic, df = degrees_freedom),
    "two.sided" = 1 - 2 * pt(-abs(statistic), df = degrees_freedom)
  )

  list(
    success = success,
    estimate = estimate,
    std_error = sqrt(total_variance),
    degrees_freedom = degrees_freedom
  )
}

#' @title Validate complete binary outcomes
#'
#' @description Ensures a binary endpoint analysis is supplied only event
#'   indicators and that censored subjects have complete follow-up.
#'
#' @noRd
assert_complete_binary_outcomes <- function(data, end_of_study, method_label) {
  tryCatch(
    validate_binary_vector(data$event, "data$event"),
    error = function(error) {
      stop(
        method_label,
        " analysis requires binary event outcomes: ",
        conditionMessage(error)
      )
    }
  )
  validate_nonnegative_numeric_vector(data$time, "data$time")
  if (any(data$event == 0 & data$time < end_of_study)) {
    stop(
      method_label,
      " analysis requires all censored subjects to be followed ",
      "to 'end_of_study' or imputed before analysis"
    )
  }
}
