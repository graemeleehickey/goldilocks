#' @title Estimate a frequentist binary risk difference
#'
#' @description Calculates the treatment-minus-control difference in event
#'   proportions and its unpooled binomial variance.
#'
#' @inheritParams analyse_data
#'
#' @return A list containing the risk-difference estimate and variance.
#'
#' @noRd
risk_difference_estimate_checked <- function(data, end_of_study) {
  assert_complete_binary_outcomes(data, end_of_study, "Risk-difference")

  treatment_event <- data$event[data$treatment == 1]
  control_event <- data$event[data$treatment == 0]
  if (length(treatment_event) == 0 || length(control_event) == 0) {
    stop("Risk-difference analysis requires at least one subject per arm")
  }

  risk_difference_from_counts(
    events_control = sum(control_event),
    n_control = length(control_event),
    events_treatment = sum(treatment_event),
    n_treatment = length(treatment_event)
  )
}

#' Estimate a binary risk difference from arm totals
#'
#' @param events_control A non-negative integer giving the number of control-arm
#'   events.
#' @param n_control A positive integer giving the number of control-arm
#'   subjects.
#' @param events_treatment A non-negative integer giving the number of
#'   treatment-arm events.
#' @param n_treatment A positive integer giving the number of treatment-arm
#'   subjects.
#'
#' @return A list containing the treatment-minus-control event-risk difference
#'   (`estimate`) and its unpooled binomial variance (`variance`).
#'
#' @keywords internal
#' @noRd
risk_difference_from_counts <- function(
  events_control,
  n_control,
  events_treatment,
  n_treatment
) {
  treatment_risk <- events_treatment / n_treatment
  control_risk <- events_control / n_control
  estimate <- treatment_risk - control_risk
  variance <- treatment_risk *
    (1 - treatment_risk) /
    n_treatment +
    control_risk * (1 - control_risk) / n_control

  if (!is.finite(estimate) || !is.finite(variance) || variance < 0) {
    stop(
      "Risk-difference analysis is non-estimable: the estimate or variance is ",
      "not finite and non-negative."
    )
  }

  list(estimate = estimate, variance = variance)
}

#' @title Test a frequentist binary risk difference
#'
#' @description Performs a Wald test of the treatment-minus-control event-risk
#'   difference against the supplied null value.
#'
#' @inheritParams analyse_data
#'
#' @return A list containing the success score, risk-difference estimate,
#'   variance, and standard error.
#'
#' @noRd
risk_difference_wald_test_checked <- function(
  data,
  end_of_study,
  alternative,
  h0
) {
  fit <- risk_difference_estimate_checked(data, end_of_study)
  risk_difference_wald_from_estimate(fit, alternative, h0)
}

#' Test a binary risk difference from arm totals
#'
#' @inheritParams risk_difference_from_counts
#' @param alternative A character value giving the direction of the alternative
#'   hypothesis: `"less"`, `"greater"`, or `"two.sided"`.
#' @param h0 A finite numeric value giving the null treatment-minus-control
#'   event-risk difference.
#'
#' @return A list containing the success score (`success`), estimated risk
#'   difference (`estimate`), variance, and standard error.
#'
#' @keywords internal
#' @noRd
risk_difference_wald_from_counts <- function(
  events_control,
  n_control,
  events_treatment,
  n_treatment,
  alternative,
  h0
) {
  fit <- risk_difference_from_counts(
    events_control = events_control,
    n_control = n_control,
    events_treatment = events_treatment,
    n_treatment = n_treatment
  )
  risk_difference_wald_from_estimate(fit, alternative, h0)
}

#' Test a binary risk difference from its estimate and variance
#'
#' @param fit A list containing a finite risk-difference estimate and a
#'   non-negative variance.
#' @inheritParams risk_difference_wald_from_counts
#'
#' @return A list containing the success score (`success`), estimated risk
#'   difference (`estimate`), variance, and standard error.
#'
#' @keywords internal
#' @noRd
risk_difference_wald_from_estimate <- function(fit, alternative, h0) {
  if (fit$variance == 0) {
    stop(
      "Risk-difference analysis is non-estimable: the estimated variance is ",
      "zero because both arms have no outcome variation."
    )
  }

  std_error <- sqrt(fit$variance)
  statistic <- (fit$estimate - h0) / std_error
  success <- normal_test_success(statistic, alternative)

  list(
    success = success,
    estimate = fit$estimate,
    variance = fit$variance,
    std_error = std_error
  )
}

#' Test a risk difference using the Farrington-Manning score statistic
#'
#' @description Validates patient-level binary outcomes and tests the
#'   treatment-minus-control risk difference using the binomial risks estimated
#'   under the null constraint.
#'
#' @inheritParams risk_difference_wald_test_checked
#'
#' @return A list containing the success score, risk-difference estimate,
#'   null-constrained variance, standard error, and score statistic.
#'
#' @noRd
risk_difference_fm_test_checked <- function(
  data,
  end_of_study,
  alternative,
  h0
) {
  assert_complete_binary_outcomes(data, end_of_study, "Risk-difference")

  treatment_event <- data$event[data$treatment == 1]
  control_event <- data$event[data$treatment == 0]
  if (length(treatment_event) == 0 || length(control_event) == 0) {
    stop("Risk-difference analysis requires at least one subject per arm")
  }

  risk_difference_fm_from_counts(
    events_control = sum(control_event),
    n_control = length(control_event),
    events_treatment = sum(treatment_event),
    n_treatment = length(treatment_event),
    alternative = alternative,
    h0 = h0
  )
}

#' Test a risk difference from arm totals using Farrington-Manning
#'
#' @description Calculates the Farrington-Manning score statistic using maximum
#'   likelihood estimates of the two event risks under the constraint that the
#'   treatment-minus-control risk difference equals `h0`.
#'
#' @inheritParams risk_difference_wald_from_counts
#'
#' @return A list containing the success score (`success`), observed risk
#'   difference (`estimate`), null-constrained variance, standard error, score
#'   statistic, and null-constrained arm risks.
#'
#' @keywords internal
#' @noRd
risk_difference_fm_from_counts <- function(
  events_control,
  n_control,
  events_treatment,
  n_treatment,
  alternative,
  h0
) {
  observed <- risk_difference_from_counts(
    events_control = events_control,
    n_control = n_control,
    events_treatment = events_treatment,
    n_treatment = n_treatment
  )
  null_risks <- risk_difference_null_risks(
    events_control = events_control,
    n_control = n_control,
    events_treatment = events_treatment,
    n_treatment = n_treatment,
    h0 = h0
  )
  variance <-
    null_risks$treatment *
    (1 - null_risks$treatment) /
    n_treatment +
    null_risks$control * (1 - null_risks$control) / n_control

  difference_from_null <- observed$estimate - h0
  statistic <- if (variance == 0) {
    if (difference_from_null == 0) 0 else sign(difference_from_null) * Inf
  } else {
    difference_from_null / sqrt(variance)
  }

  list(
    success = normal_test_success(statistic, alternative),
    estimate = observed$estimate,
    variance = variance,
    std_error = sqrt(variance),
    statistic = statistic,
    null_control_risk = null_risks$control,
    null_treatment_risk = null_risks$treatment
  )
}

#' Estimate event risks under a null risk-difference constraint
#'
#' @description Maximizes the joint binomial likelihood subject to the treatment
#'   risk minus the control risk equaling `h0`.
#'
#' @inheritParams risk_difference_fm_from_counts
#'
#' @return A list containing the null-constrained control and treatment risks.
#'
#' @keywords internal
#' @noRd
risk_difference_null_risks <- function(
  events_control,
  n_control,
  events_treatment,
  n_treatment,
  h0
) {
  lower_control <- max(0, -h0)
  upper_control <- min(1, 1 - h0)

  if (h0 == 0) {
    control_risk <-
      (events_control + events_treatment) / (n_control + n_treatment)
  } else if (lower_control == upper_control) {
    control_risk <- lower_control
  } else {
    log_likelihood <- function(control_risk) {
      treatment_risk <- control_risk + h0
      stats::dbinom(
        events_control,
        size = n_control,
        prob = control_risk,
        log = TRUE
      ) +
        stats::dbinom(
          events_treatment,
          size = n_treatment,
          prob = treatment_risk,
          log = TRUE
        )
    }
    optimum <- stats::optimize(
      function(control_risk) -log_likelihood(control_risk),
      interval = c(lower_control, upper_control),
      tol = .Machine$double.eps^0.75
    )$minimum
    candidates <- c(lower_control, optimum, upper_control)
    likelihoods <- vapply(candidates, log_likelihood, numeric(1))
    control_risk <- candidates[[which.max(likelihoods)]]
  }

  list(
    control = control_risk,
    treatment = control_risk + h0
  )
}
