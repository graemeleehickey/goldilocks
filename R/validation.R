validate_single_probability <- function(x, name, upper_open = FALSE) {
  if (
    length(x) != 1 ||
      !is.numeric(x) ||
      is.na(x) ||
      !is.finite(x) ||
      x < 0 ||
      x > 1 ||
      (upper_open && x >= 1)
  ) {
    bound <- if (upper_open) "[0, 1)" else "[0, 1]"
    stop("'", name, "' must be a single finite probability in ", bound)
  }

  invisible(TRUE)
}

validate_probability_vector <- function(x, name, upper_open = FALSE) {
  if (
    !is.numeric(x) ||
      length(x) == 0 ||
      any(is.na(x)) ||
      any(!is.finite(x)) ||
      any(x < 0) ||
      any(x > 1) ||
      (upper_open && any(x >= 1))
  ) {
    bound <- if (upper_open) "[0, 1)" else "[0, 1]"
    stop("'", name, "' must contain finite probabilities in ", bound)
  }

  invisible(TRUE)
}

validate_positive_integer_scalar <- function(x, name) {
  if (
    length(x) != 1 ||
      !is.numeric(x) ||
      is.na(x) ||
      !is.finite(x) ||
      x <= 0 ||
      x != floor(x)
  ) {
    stop("'", name, "' must be a single positive integer")
  }

  invisible(TRUE)
}

validate_positive_integer_vector <- function(x, name) {
  if (
    !is.numeric(x) ||
      length(x) == 0 ||
      any(is.na(x)) ||
      any(!is.finite(x)) ||
      any(x <= 0) ||
      any(x != floor(x))
  ) {
    stop("'", name, "' must contain positive integers")
  }

  invisible(TRUE)
}

validate_gamma_prior <- function(prior, name = "prior") {
  if (
    !is.numeric(prior) ||
      length(prior) != 2 ||
      any(is.na(prior)) ||
      any(!is.finite(prior)) ||
      any(prior <= 0)
  ) {
    stop("'", name, "' must contain two positive finite values")
  }

  invisible(TRUE)
}

validate_cutpoints <- function(cutpoints) {
  if (
    !is.numeric(cutpoints) ||
      length(cutpoints) == 0 ||
      any(is.na(cutpoints)) ||
      any(!is.finite(cutpoints))
  ) {
    stop("'cutpoints' must contain finite numeric values")
  }

  if (cutpoints[1] != 0) {
    stop("First element of 'cutpoints' should be 0")
  }

  # A non-increasing knot vector creates zero-width or backward intervals,
  # which makes piecewise exposure and event allocation undefined.
  if (any(diff(cutpoints) <= 0)) {
    stop("'cutpoints' must be strictly increasing")
  }

  invisible(TRUE)
}

validate_endpoint_time <- function(endpoint, cutpoints, name, after_last = TRUE) {
  if (
    length(endpoint) != 1 ||
      !is.numeric(endpoint) ||
      is.na(endpoint) ||
      !is.finite(endpoint) ||
      endpoint <= 0
  ) {
    stop("'", name, "' must be a single finite positive value")
  }

  if (after_last && endpoint <= max(cutpoints)) {
    stop("'", name, "' must be a finite value greater than the last cutpoint")
  }

  invisible(TRUE)
}

validate_piecewise_hazard <- function(hazard, cutpoints, name = "hazard") {
  if (
    !is.numeric(hazard) ||
      !is.null(dim(hazard)) ||
      length(hazard) == 0 ||
      any(is.na(hazard)) ||
      any(!is.finite(hazard)) ||
      any(hazard < 0)
  ) {
    stop("'", name, "' must contain finite non-negative hazard rates")
  }

  if (length(hazard) != length(cutpoints)) {
    stop("Length of 'cutpoints' must be equal to length of '", name, "'")
  }

  invisible(TRUE)
}

validate_hazard_matrix <- function(hazard, cutpoints, name = "hazard") {
  if (
    !is.matrix(hazard) ||
      !is.numeric(hazard) ||
      nrow(hazard) == 0 ||
      any(is.na(hazard)) ||
      any(!is.finite(hazard)) ||
      any(hazard < 0)
  ) {
    stop("'", name, "' must be a non-empty matrix of finite non-negative hazard rates")
  }

  if (ncol(hazard) != length(cutpoints)) {
    stop("The length of the hazard rates and cutpoints do not match")
  }

  invisible(TRUE)
}

validate_enrollment_schedule <- function(lambda, lambda_time, N_total) {
  validate_positive_integer_scalar(N_total, "N_total")

  # A non-positive or non-finite rate can leave the accrual loop unable to
  # reach its target sample size, so reject it before any random draws occur.
  if (
    !is.numeric(lambda) ||
      !is.null(dim(lambda)) ||
      length(lambda) == 0 ||
      any(is.na(lambda)) ||
      any(!is.finite(lambda)) ||
      any(lambda <= 0)
  ) {
    stop("'lambda' must contain finite positive enrollment rates")
  }

  if (length(lambda) != length(lambda_time)) {
    stop("The length of rates should match the length of knots")
  }

  if (
    !is.numeric(lambda_time) ||
      !is.null(dim(lambda_time)) ||
      length(lambda_time) == 0 ||
      any(is.na(lambda_time)) ||
      any(!is.finite(lambda_time))
  ) {
    stop("'lambda_time' must contain finite numeric knot values")
  }

  if (lambda_time[1] != 0) {
    stop("The first cutpoint should always 0")
  }

  # The rate-selection loop assumes chronological knots. Without this check,
  # an unordered schedule can silently assign subjects to the wrong rate.
  if (any(diff(lambda_time) <= 0)) {
    stop("'lambda_time' must be strictly increasing")
  }

  invisible(TRUE)
}

validate_h0 <- function(h0, method, single_arm) {
  if (
    length(h0) != 1 ||
      !is.numeric(h0) ||
      is.na(h0) ||
      !is.finite(h0)
  ) {
    stop("'h0' must be a single finite numeric value")
  }

  if (method %in% c("bayes-surv", "bayes-bin")) {
    lower <- if (single_arm) 0 else -1
    upper <- 1
    if (h0 < lower || h0 > upper) {
      stop(
        "'h0' must lie in [",
        lower,
        ", ",
        upper,
        "] for ",
        if (single_arm) "single-arm" else "two-arm",
        " Bayesian analyses"
      )
    }
  }

  invisible(TRUE)
}
