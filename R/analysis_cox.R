#' Calculate a Cox Wald result from outcome vectors
#'
#' @inheritParams analyse_logrank
#' @param h0 A finite numeric value giving the null log hazard ratio.
#'
#' @return A list containing `1 - P` for the Cox Wald test and the estimated log
#'   hazard ratio.
#'
#' @keywords internal
#' @noRd
analyse_cox <- function(time, event, treatment, alternative, h0) {
  fit_cox <- cox_wald_outcomes_checked(time, event, treatment)
  z <- (fit_cox$estimate - h0) / fit_cox$std_error
  if (alternative == "two.sided") {
    success <- 1 - (2 * pnorm(-abs(z)))
  } else {
    # Cox z < 0 when the estimated log hazard ratio is less than h0.
    # "less" => treatment beneficial/non-inferior => large success when z << 0
    # "greater" => treatment harmful/superior to h0 => large success when z >> 0
    if (alternative == "less") {
      success <- 1 - pnorm(z)
    } else if (alternative == "greater") {
      success <- pnorm(z)
    }
  }

  list(success = success, effect = fit_cox$estimate)
}

#' @title Calculate a Cox proportional hazards Wald test
#'
#' @description Fits a Cox proportional hazards model for the treatment effect
#'   and converts fitting warnings into clear non-estimability errors.
#'
#' @inheritParams analyse_data
#'
#' @return A list with the estimated log hazard ratio (`estimate`) and its
#'   standard error (`std_error`).
#'
#' @noRd
cox_wald_test_checked <- function(data, ...) {
  cox_wald_outcomes_checked(
    time = data$time,
    event = data$event,
    treatment = data$treatment,
    ...
  )
}

#' Calculate a checked Cox Wald test from outcome vectors
#'
#' @inheritParams analyse_logrank
#'
#' @return A list with the estimated log hazard ratio (`estimate`) and its
#'   standard error (`std_error`).
#'
#' @keywords internal
#' @noRd
cox_wald_outcomes_checked <- function(time, event, treatment, ...) {
  fit_state <- new.env(parent = emptyenv())
  fit_state$warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(
      cox_wald_from_outcomes(time, event, treatment, ...),
      warning = function(w) {
        fit_state$warnings <- c(
          fit_state$warnings,
          conditionMessage(w)
        )
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      stop(
        "Cox analysis is non-estimable: the Cox fitter failed: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  if (length(fit_state$warnings) > 0L) {
    stop(
      "Cox analysis is non-estimable: the Cox model did not produce a ",
      "reliable Wald statistic. The Cox fitter reported: ",
      paste(unique(fit_state$warnings), collapse = " | "),
      call. = FALSE
    )
  }

  assert_cox_estimable(fit)
  fit
}

#' @title Assert that a Cox result is estimable
#'
#' @description Ensures a Cox treatment-effect estimate and its standard error
#'   are finite before they are used in an adaptive decision.
#'
#' @noRd
assert_cox_estimable <- function(fit) {
  if (
    length(fit$estimate) != 1 ||
      length(fit$std_error) != 1 ||
      is.na(fit$estimate) ||
      is.na(fit$std_error) ||
      !is.finite(fit$estimate) ||
      !is.finite(fit$std_error) ||
      fit$std_error <= 0
  ) {
    stop(
      "Cox analysis is non-estimable: the treatment effect or standard ",
      "error is not finite and positive. This can occur when there are no ",
      "events, separation, or insufficient information to estimate the ",
      "treatment effect."
    )
  }

  invisible(TRUE)
}

#' @title Fit a Cox Wald test using an available survival calculation
#'
#' @description Uses the lower-overhead `coxph.fit()` calculation when its
#'   arguments are compatible with the installed version of `survival`, and
#'   otherwise uses [survival::coxph()]. The `engine` and `compatibility`
#'   arguments support equivalence testing of the two calculations.
#'
#' @noRd
cox_wald_test <- function(
  data,
  engine = c("auto", "fast", "public"),
  compatibility = coxph_fit_compatibility()
) {
  cox_wald_from_outcomes(
    time = data$time,
    event = data$event,
    treatment = data$treatment,
    engine = engine,
    compatibility = compatibility
  )
}

#' Fit a Cox Wald test from outcome vectors
#'
#' @inheritParams analyse_logrank
#' @param engine A character value selecting the automatic, faster, or exported
#'   Cox calculation.
#' @param compatibility The result of `coxph_fit_compatibility()`.
#'
#' @return A list with the estimated log hazard ratio (`estimate`) and its
#'   standard error (`std_error`).
#'
#' @keywords internal
#' @noRd
cox_wald_from_outcomes <- function(
  time,
  event,
  treatment,
  engine = c("auto", "fast", "public"),
  compatibility = coxph_fit_compatibility()
) {
  engine <- match.arg(engine)
  if (engine == "public") {
    return(cox_wald_with_coxph(time, event, treatment))
  }

  fast_available <- isTRUE(compatibility$compatible) &&
    is.function(compatibility$fitter)
  if (!fast_available) {
    if (engine == "fast") {
      stop(
        "The survival::coxph.fit() fast path is unavailable: ",
        compatibility$reason,
        call. = FALSE
      )
    }
    return(cox_wald_with_coxph(time, event, treatment))
  }

  fast_fit <- tryCatch(
    cox_wald_with_coxph_fit(time, event, treatment, compatibility$fitter),
    error = identity
  )
  if (!inherits(fast_fit, "error")) {
    return(fast_fit)
  }
  if (engine == "fast") {
    stop(fast_fit)
  }

  cox_wald_with_coxph(time, event, treatment)
}

#' @title Check availability of the lower-overhead Cox calculation
#'
#' @description Records the installed `survival` version and verifies every
#'   unexported fitter argument used by the package. Installations outside the
#'   supported compatibility boundary use the exported fallback.
#'
#' @noRd
coxph_fit_compatibility <- local({
  cached <- NULL

  function(refresh = FALSE) {
    if (!refresh && !is.null(cached)) {
      return(cached)
    }

    version <- as.character(utils::packageVersion("survival"))
    fitter <- tryCatch(
      getFromNamespace("coxph.fit", "survival"),
      error = function(e) NULL
    )
    required_arguments <- c(
      "x",
      "y",
      "strata",
      "offset",
      "init",
      "control",
      "weights",
      "method",
      "rownames",
      "resid",
      "nocenter"
    )
    version_supported <- utils::compareVersion(version, "3.2-0") >= 0L
    missing_arguments <- if (is.function(fitter)) {
      setdiff(required_arguments, names(formals(fitter)))
    } else {
      required_arguments
    }
    compatible <- version_supported &&
      is.function(fitter) &&
      length(missing_arguments) == 0L
    reason <- if (!version_supported) {
      paste0("survival ", version, " predates the guarded fast path")
    } else if (!is.function(fitter)) {
      "survival::coxph.fit() could not be resolved"
    } else if (length(missing_arguments) > 0L) {
      paste0(
        "survival::coxph.fit() lacks required argument(s): ",
        paste(missing_arguments, collapse = ", ")
      )
    } else {
      paste0("compatible survival ", version, " signature")
    }

    cached <<- list(
      compatible = compatible,
      fitter = if (compatible) fitter else NULL,
      version = version,
      reason = reason
    )
    cached
  }
})

#' Fit a Cox Wald test from outcome vectors with `coxph.fit()`
#'
#' @inheritParams analyse_logrank
#' @param fitter A compatible `survival::coxph.fit()` function.
#'
#' @return A list with the estimated log hazard ratio (`estimate`) and its
#'   standard error (`std_error`).
#'
#' @keywords internal
#' @noRd
cox_wald_with_coxph_fit <- function(time, event, treatment, fitter) {
  y <- survival::Surv(time, event)
  x <- matrix(
    as.double(treatment),
    ncol = 1L,
    dimnames = list(NULL, "treatment")
  )
  fit <- fitter(
    x = x,
    y = y,
    strata = NULL,
    offset = NULL,
    init = NULL,
    control = survival::coxph.control(),
    weights = NULL,
    method = "efron",
    rownames = NULL,
    resid = FALSE,
    nocenter = NULL
  )

  valid_result <- is.list(fit) &&
    is.numeric(fit$coefficients) &&
    identical(names(fit$coefficients), "treatment") &&
    is.matrix(fit$var) &&
    identical(dim(fit$var), c(1L, 1L))
  if (!valid_result) {
    stop(
      "survival::coxph.fit() returned an incompatible result structure",
      call. = FALSE
    )
  }

  list(
    estimate = unname(fit$coefficients[["treatment"]]),
    std_error = sqrt(unname(fit$var[1L, 1L]))
  )
}

#' Fit a Cox Wald test from outcome vectors with [survival::coxph()]
#'
#' @inheritParams analyse_logrank
#'
#' @return A list with the estimated log hazard ratio (`estimate`) and its
#'   standard error (`std_error`).
#'
#' @keywords internal
#' @noRd
cox_wald_with_coxph <- function(time, event, treatment) {
  data <- data.frame(
    time = time,
    event = event,
    treatment = treatment
  )
  fit <- survival::coxph(
    survival::Surv(time, event) ~ treatment,
    data = data,
    ties = "efron",
    singular.ok = FALSE,
    model = FALSE,
    x = FALSE,
    y = FALSE
  )
  estimate <- stats::coef(fit)["treatment"]
  variance <- stats::vcov(fit)["treatment", "treatment"]

  list(
    estimate = unname(estimate),
    std_error = sqrt(unname(variance))
  )
}
