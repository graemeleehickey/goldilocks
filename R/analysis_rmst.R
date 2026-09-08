#' Validate a prespecified RMST horizon and null difference
#'
#' @inheritParams survival_adapt
#' @noRd
validate_rmst_args <- function(rmst_tau, end_of_study, h0) {
  validate_endpoint_time(rmst_tau, NULL, "rmst_tau")
  if (rmst_tau > end_of_study) {
    stop("'rmst_tau' must not exceed 'end_of_study'", call. = FALSE)
  }
  validate_h0(h0, "rmst", FALSE)
  if (abs(h0) > rmst_tau) {
    stop(
      "'h0' must lie in [-rmst_tau, rmst_tau] for RMST analyses",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Estimate a treatment-control restricted mean survival time difference
#'
#' @description Integrates each Kaplan-Meier curve from zero to the fixed
#'   horizon, using the Greenwood plug-in variance. Zero within-arm or
#'   within-imputation variances are retained for subsequent pooling; a test
#'   requires positive total variance.
#' @inheritParams analyse_logrank
#' @inheritParams survival_adapt
#' @param engine Internal calculation selector. `"auto"` uses the equivalent
#'   empirical calculation when every truncated outcome is known; `"survfit"`
#'   always uses Kaplan-Meier integration for reference checks.
#'
#' @return A list with `estimate`, `variance`, `std_error`, and named `rmst` and
#'   `arm_variance` vectors (control, treatment).
#' @noRd
rmst_estimate <- function(
  time,
  event,
  treatment,
  rmst_tau,
  engine = c("auto", "survfit")
) {
  engine <- match.arg(engine)
  validate_endpoint_time(rmst_tau, NULL, "rmst_tau")
  validate_nonnegative_numeric_vector(time, "time")
  validate_binary_vector(event, "event")
  validate_binary_vector(treatment, "treatment")
  if (length(time) != length(event) || length(time) != length(treatment)) {
    stop("RMST analysis requires equally sized outcome vectors", call. = FALSE)
  }
  if (!all(c(0, 1) %in% treatment)) {
    stop(
      "RMST analysis is non-estimable: both treatment arms are required",
      call. = FALSE
    )
  }

  if (engine == "auto" && all(event == 1 | time >= rmst_tau)) {
    # With no censoring before tau, KM integration equals mean(min(T, tau)).
    # Greenwood gives sum((x - mean(x))^2) / n^2, not var(x) / n.
    # Predictively completed outcomes always meet this condition.
    x <- split(pmin(time, rmst_tau), as.integer(treatment))
    rmst <- stats::setNames(
      vapply(x, mean, numeric(1)),
      c("control", "treatment")
    )
    arm_variance <- stats::setNames(
      vapply(
        x,
        function(v) {
          sum((v - mean(v))^2) / length(v)^2
        },
        numeric(1)
      ),
      c("control", "treatment")
    )
  } else {
    # Explicit numeric coding gives stable control/treatment ordering, including
    # when the input indicators are logical. Censoring at tau does not alter
    # the area or its variance; events exactly at tau contribute zero tail area.
    data <- data.frame(
      time = pmin(time, rmst_tau),
      event = as.integer(event == 1 & time <= rmst_tau),
      treatment = as.integer(treatment)
    )
    fit <- survival::survfit(
      survival::Surv(time, event) ~ treatment,
      data = data
    )
    arm_ends <- cumsum(fit$strata)
    followup_ends <- vapply(split(data$time, data$treatment), max, numeric(1))
    unsupported <- followup_ends < rmst_tau & fit$surv[arm_ends] > 0
    if (any(unsupported)) {
      stop(
        "RMST analysis is non-estimable: follow-up ends before 'rmst_tau' ",
        "with positive survival in the ",
        paste(c("control", "treatment")[unsupported], collapse = " and "),
        " arm; the fixed horizon cannot be estimated without extrapolation",
        call. = FALSE
      )
    }

    tab <- summary(fit, rmean = rmst_tau)$table
    rmst <- stats::setNames(
      as.numeric(tab[, "rmean"]),
      c("control", "treatment")
    )
    arm_variance <- stats::setNames(
      as.numeric(tab[, "se(rmean)"])^2,
      c("control", "treatment")
    )
  }
  if (
    any(!is.finite(rmst)) ||
      any(!is.finite(arm_variance)) ||
      any(arm_variance < 0)
  ) {
    stop(
      "RMST analysis is non-estimable: non-finite estimate or variance",
      call. = FALSE
    )
  }
  variance <- sum(arm_variance)
  list(
    estimate = unname(rmst["treatment"] - rmst["control"]),
    variance = variance,
    std_error = sqrt(variance),
    rmst = rmst,
    arm_variance = arm_variance
  )
}

#' Require positive information for an RMST Wald test
#'
#' @param variance The variance of an RMST difference, possibly Rubin-pooled.
#' @noRd
assert_rmst_variance <- function(variance) {
  if (length(variance) != 1L || !is.finite(variance) || variance <= 0) {
    stop(
      "RMST analysis is non-estimable: total variance must be positive",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Apply an RMST difference Wald test to completed outcomes
#'
#' @inheritParams rmst_estimate
#' @inheritParams survival_adapt
#' @return The completed-data `success` and `effect` contract.
#' @noRd
analyse_rmst <- function(time, event, treatment, rmst_tau, alternative, h0) {
  validate_rmst_args(rmst_tau, rmst_tau, h0)
  fit <- rmst_estimate(time, event, treatment, rmst_tau)
  assert_rmst_variance(fit$variance)
  list(
    success = normal_test_success(
      (fit$estimate - h0) / fit$std_error,
      alternative
    ),
    effect = fit$estimate
  )
}
