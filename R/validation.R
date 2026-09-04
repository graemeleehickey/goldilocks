#' @title Validate one probability
#'
#' @description Checks that an input is a single finite probability, optionally
#'   excluding one to support distributions with an open upper bound.
#'
#' @param x A numeric value to validate.
#' @param name A single character string naming the argument in error messages.
#' @param upper_open A single logical value indicating whether `1` should be
#'   excluded. The default is `FALSE`.
#'
#' @noRd
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

#' @title Validate a probability vector
#'
#' @description Checks that every value in an input vector is a finite
#'   probability, optionally excluding one.
#'
#' @param x A numeric vector to validate.
#' @param name A single character string naming the argument in error messages.
#' @param upper_open A single logical value indicating whether `1` should be
#'   excluded. The default is `FALSE`.
#'
#' @noRd
validate_probability_vector <- function(x, name, upper_open = FALSE) {
  if (
    !is.numeric(x) ||
      !is.null(dim(x)) ||
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

#' Normalize a scalar or two-arm vector
#'
#' @description Converts supported arm-vector inputs into the canonical order
#'   `control`, `treatment`. Callers select whether a shared scalar or an
#'   unnamed legacy vector is allowed.
#'
#' @param x A vector containing a shared or arm-specific value.
#' @param name A single character string naming the argument in error messages.
#' @param single_arm A single logical value indicating whether the design has
#'   only a treatment arm. The default is `FALSE`.
#' @param allow_scalar A single logical value indicating whether a shared value
#'   may be repeated across two arms. The default is `FALSE`.
#' @param allow_unnamed A single logical value indicating whether an unnamed
#'   length-two vector may be interpreted in legacy `c(control, treatment)`
#'   order. The default is `FALSE`.
#' @param warn_unnamed_unequal A single logical value indicating whether an
#'   unequal unnamed vector should produce a future-compatibility warning. The
#'   default is `FALSE`.
#'
#' @return A named vector containing `treatment` for a single-arm design or
#'   `control` followed by `treatment` for a two-arm design.
#'
#' @keywords internal
#' @noRd
normalize_arm_vector <- function(
  x,
  name,
  single_arm = FALSE,
  allow_scalar = FALSE,
  allow_unnamed = FALSE,
  warn_unnamed_unequal = FALSE
) {
  if (single_arm) {
    if (length(x) != 1L) {
      stop("'", name, "' must be a single value for a single-arm design")
    }
    return(c(treatment = unname(x)))
  }

  arm_names <- c("control", "treatment")
  if (length(x) == 1L && allow_scalar) {
    normalized <- rep(unname(x), 2L)
    names(normalized) <- arm_names
    return(normalized)
  }

  supplied_names <- names(x)
  valid_named_vector <- length(x) == 2L &&
    !is.null(supplied_names) &&
    !anyNA(supplied_names) &&
    all(nzchar(supplied_names)) &&
    !anyDuplicated(supplied_names) &&
    setequal(supplied_names, arm_names)
  if (valid_named_vector) {
    normalized <- unname(x[arm_names])
    names(normalized) <- arm_names
    return(normalized)
  }

  fully_unnamed <- is.null(supplied_names) ||
    (!anyNA(supplied_names) && all(!nzchar(supplied_names)))
  if (length(x) == 2L && allow_unnamed && fully_unnamed) {
    if (warn_unnamed_unequal && length(unique(unname(x))) > 1L) {
      warning(
        "Using unequal unnamed '",
        name,
        "' assumes c(control, treatment) order; name the values explicitly ",
        "because unequal unnamed arm allocations may require names in a ",
        "future major release",
        call. = FALSE
      )
    }
    normalized <- unname(x)
    names(normalized) <- arm_names
    return(normalized)
  }

  scalar_text <- if (allow_scalar) "a single value or " else ""
  unnamed_text <- if (allow_unnamed) {
    "an unnamed length-two vector in c(control, treatment) order or "
  } else {
    ""
  }
  stop(
    "'",
    name,
    "' must be ",
    scalar_text,
    unnamed_text,
    "a length-two vector named 'control' and 'treatment' for a two-arm design"
  )
}

#' Normalize arm-specific loss-to-follow-up proportions
#'
#' @param prop_loss A numeric vector containing one shared probability or two
#'   values named `control` and `treatment`.
#' @param single_arm A single logical value indicating whether the design has
#'   only a treatment arm.
#'
#' @return A canonical named arm vector.
#'
#' @keywords internal
#' @noRd
normalize_prop_loss <- function(prop_loss, single_arm) {
  validate_probability_vector(prop_loss, "prop_loss")
  if (single_arm && length(prop_loss) != 1L) {
    stop("'prop_loss' must be a single probability for a single-arm design")
  }
  normalize_arm_vector(
    prop_loss,
    name = "prop_loss",
    single_arm = single_arm,
    allow_scalar = TRUE
  )
}

#' Normalize a two-arm randomization allocation
#'
#' @param allocation A length-two positive integer vector giving the control and
#'   treatment allocation weights.
#' @param name A single character string naming the argument in error messages.
#'   The default is `"allocation"`.
#'
#' @return A positive integer vector in canonical arm order.
#'
#' @keywords internal
#' @noRd
normalize_allocation <- function(allocation, name = "allocation") {
  if (length(allocation) != 2L) {
    stop("'", name, "' must contain two positive integer values")
  }
  if (
    !is.numeric(allocation) ||
      !is.null(dim(allocation)) ||
      any(is.na(allocation)) ||
      any(!is.finite(allocation)) ||
      any(allocation %% 1 != 0)
  ) {
    stop("All values of '", name, "' must be integer values")
  }
  if (any(allocation <= 0)) {
    stop("'", name, "' must contain two positive integer values")
  }

  normalize_arm_vector(
    allocation,
    name = name,
    allow_unnamed = TRUE,
    warn_unnamed_unequal = TRUE
  )
}

#' @title Validate one positive integer
#'
#' @description Checks that an input is a single finite, strictly positive
#'   integer suitable for a sample size or count.
#'
#' @param x A numeric value to validate.
#' @param name A single character string naming the argument in error messages.
#'
#' @noRd
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

#' @title Validate one non-negative integer
#'
#' @description Checks that an input is a single finite, non-negative integer.
#'   Unlike `validate_positive_integer_scalar()`, zero is accepted so callers
#'   can define an explicit zero-draw result.
#'
#' @param x A numeric value to validate.
#' @param name A single character string naming the argument in error messages.
#'
#' @noRd
validate_nonnegative_integer_scalar <- function(x, name) {
  if (
    length(x) != 1 ||
      !is.numeric(x) ||
      is.na(x) ||
      !is.finite(x) ||
      x < 0 ||
      x != floor(x)
  ) {
    stop("'", name, "' must be a single non-negative integer")
  }

  invisible(TRUE)
}

#' @title Validate a positive integer vector
#'
#' @description Checks that an input vector contains only finite, strictly
#'   positive integers.
#'
#' @param x A numeric vector to validate.
#' @param name A single character string naming the argument in error messages.
#'
#' @noRd
validate_positive_integer_vector <- function(x, name) {
  if (
    !is.numeric(x) ||
      !is.null(dim(x)) ||
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

#' @title Validate one logical value
#'
#' @description Checks that an input is a non-missing scalar logical value.
#'
#' @param x A logical value to validate.
#' @param name A single character string naming the argument in error messages.
#'
#' @noRd
validate_logical_scalar <- function(x, name) {
  if (!is.logical(x) || length(x) != 1 || is.na(x)) {
    stop("'", name, "' must be TRUE or FALSE")
  }

  invisible(TRUE)
}

#' @title Validate finite non-negative values
#'
#' @description Checks that a numeric vector contains finite, non-negative
#'   values. Callers opt in explicitly when a zero-length vector has a defined
#'   result.
#'
#' @param x A numeric vector to validate.
#' @param name A single character string naming the argument in error messages.
#' @param allow_empty A single logical value indicating whether a zero-length
#'   vector is permitted. The default is `FALSE`.
#'
#' @noRd
validate_nonnegative_numeric_vector <- function(
  x,
  name,
  allow_empty = FALSE
) {
  if (
    !is.numeric(x) ||
      !is.null(dim(x)) ||
      (!allow_empty && length(x) == 0) ||
      anyNA(x) ||
      any(!is.finite(x)) ||
      any(x < 0)
  ) {
    empty_note <- if (allow_empty) "" else " non-empty"
    stop(
      "'",
      name,
      "' must be a",
      empty_note,
      " numeric vector of finite non-negative values"
    )
  }

  invisible(TRUE)
}

#' @title Validate binary indicators
#'
#' @description Checks that a vector contains only non-missing zero-one event
#'   indicators.
#'
#' @param x A numeric or logical vector to validate.
#' @param name A single character string naming the argument in error messages.
#' @param allow_empty A single logical value indicating whether a zero-length
#'   vector is permitted. The default is `FALSE`.
#'
#' @noRd
validate_binary_vector <- function(x, name, allow_empty = FALSE) {
  if (
    !(is.numeric(x) || is.logical(x)) ||
      !is.null(dim(x)) ||
      (!allow_empty && length(x) == 0) ||
      anyNA(x) ||
      any(!is.finite(x)) ||
      any(!x %in% c(0, 1))
  ) {
    empty_note <- if (allow_empty) "" else " non-empty"
    stop(
      "'",
      name,
      "' must be a",
      empty_note,
      " numeric or logical vector containing only 0 and 1"
    )
  }

  invisible(TRUE)
}

#' @title Validate and normalize a Gamma survival prior
#'
#' @description Checks that Gamma priors supply finite, strictly positive shape
#'   and rate parameters, then broadcasts shared or arm-specific inputs to a
#'   named parameter-by-interval-by-arm array.
#'
#' @param prior_surv A numeric vector, matrix, list, or three-dimensional array
#'   containing Gamma shape and rate parameters.
#' @param n_intervals A positive integer giving the number of
#'   piecewise-exponential intervals.
#' @param single_arm A single logical value indicating whether the design has
#'   only a treatment arm.
#' @param name A single character string naming the prior argument in error
#'   messages. The default is `"prior_surv"`.
#'
#' @noRd
normalize_gamma_prior <- function(
  prior_surv,
  n_intervals,
  single_arm,
  name = "prior_surv"
) {
  validate_positive_integer_scalar(n_intervals, "n_intervals")
  validate_logical_scalar(single_arm, "single_arm")
  arm_names <- if (single_arm) "treatment" else c("control", "treatment")
  interval_names <- as.character(seq_len(n_intervals))

  normalize_one_arm <- function(x, component_name) {
    valid_values <- is.numeric(x) &&
      !anyNA(x) &&
      all(is.finite(x)) &&
      all(x > 0)
    is_vector <- is.null(dim(x)) && length(x) == 2L
    is_interval_matrix <- is.matrix(x) &&
      identical(dim(x), c(2L, as.integer(n_intervals)))

    if (!valid_values || (!is_vector && !is_interval_matrix)) {
      stop(
        "'",
        component_name,
        "' must be a length-two positive finite numeric vector or a 2 x ",
        n_intervals,
        " matrix with shape in row 1 and rate in row 2"
      )
    }

    if (is_vector) {
      x <- matrix(rep(x, n_intervals), nrow = 2L)
    }
    dimnames(x) <- list(c("shape", "rate"), interval_names)
    x
  }

  # Accept the canonical representation so normalization is idempotent when a
  # validated prior passes through more than one calculation layer.
  is_prior_array <- is.array(prior_surv) &&
    length(dim(prior_surv)) == 3L
  if (is_prior_array) {
    valid_values <- is.numeric(prior_surv) &&
      !anyNA(prior_surv) &&
      all(is.finite(prior_surv)) &&
      all(prior_surv > 0)
    supplied_arms <- dimnames(prior_surv)[[3L]]
    valid_dimensions <- identical(
      dim(prior_surv),
      c(2L, as.integer(n_intervals), length(arm_names))
    )
    valid_arms <- !is.null(supplied_arms) &&
      !anyNA(supplied_arms) &&
      !anyDuplicated(supplied_arms) &&
      setequal(supplied_arms, arm_names)
    if (!valid_values || !valid_dimensions || !valid_arms) {
      stop(
        "'",
        name,
        "' array must have dimensions 2 x ",
        n_intervals,
        " x ",
        length(arm_names),
        " with finite positive values and arm names: ",
        paste(arm_names, collapse = ", ")
      )
    }
    prior_surv <- prior_surv[,, arm_names, drop = FALSE]
    dimnames(prior_surv) <- list(
      parameter = c("shape", "rate"),
      interval = interval_names,
      arm = arm_names
    )
    return(prior_surv)
  }

  if (is.list(prior_surv) && !is.data.frame(prior_surv)) {
    supplied_arms <- names(prior_surv)
    valid_arms <- !is.null(supplied_arms) &&
      !anyNA(supplied_arms) &&
      all(nzchar(supplied_arms)) &&
      !anyDuplicated(supplied_arms) &&
      setequal(supplied_arms, arm_names)
    if (!valid_arms) {
      stop(
        "'",
        name,
        "' arm-specific list must contain exactly: ",
        paste(arm_names, collapse = ", ")
      )
    }
    arm_priors <- lapply(arm_names, function(arm) {
      normalize_one_arm(prior_surv[[arm]], paste0(name, "$", arm))
    })
    names(arm_priors) <- arm_names
  } else {
    shared_prior <- normalize_one_arm(prior_surv, name)
    arm_priors <- rep(list(shared_prior), length(arm_names))
    names(arm_priors) <- arm_names
  }

  normalized <- array(
    NA_real_,
    dim = c(2L, n_intervals, length(arm_names)),
    dimnames = list(
      parameter = c("shape", "rate"),
      interval = interval_names,
      arm = arm_names
    )
  )
  for (arm in arm_names) {
    normalized[,, arm] <- arm_priors[[arm]]
  }
  normalized
}

#' @title Validate piecewise cutpoints
#'
#' @description Checks that a piecewise model has finite, positive, strictly
#'   increasing interior cutpoints, or no cutpoints for a constant hazard.
#'
#' @param cutpoints `NULL`, or a numeric vector of interior cutpoints.
#' @param name A single character string naming the argument in error messages.
#'   The default is `"cutpoints"`.
#'
#' @noRd
validate_cutpoints <- function(cutpoints, name = "cutpoints") {
  if (is.null(cutpoints)) {
    return(invisible(TRUE))
  }

  if (
    !is.numeric(cutpoints) ||
      !is.null(dim(cutpoints)) ||
      length(cutpoints) == 0 ||
      any(is.na(cutpoints)) ||
      any(!is.finite(cutpoints)) ||
      any(cutpoints <= 0)
  ) {
    stop("'", name, "' must be NULL or contain finite positive numeric values")
  }

  # A non-increasing knot vector creates zero-width or backward intervals,
  # which makes piecewise exposure and event allocation undefined.
  if (any(diff(cutpoints) <= 0)) {
    stop("'", name, "' must be strictly increasing")
  }

  invisible(TRUE)
}

#' @title Validate an endpoint time
#'
#' @description Checks that an analysis endpoint is finite and positive and
#'   lies after the final piecewise cutpoint.
#'
#' @param endpoint A numeric value giving the endpoint time.
#' @param cutpoints `NULL`, or a numeric vector of interior cutpoints.
#' @param name A single character string naming the endpoint argument in error
#'   messages.
#' @param cutpoints_name A single character string naming the cutpoint argument
#'   in error messages. The default is `"cutpoints"`.
#'
#' @noRd
validate_endpoint_time <- function(
  endpoint,
  cutpoints,
  name,
  cutpoints_name = "cutpoints"
) {
  if (
    length(endpoint) != 1 ||
      !is.numeric(endpoint) ||
      is.na(endpoint) ||
      !is.finite(endpoint) ||
      endpoint <= 0
  ) {
    stop("'", name, "' must be a single finite positive value")
  }

  if (length(cutpoints) > 0 && endpoint <= max(cutpoints)) {
    if (identical(cutpoints_name, "cutpoints")) {
      stop("'", name, "' must be a finite value greater than the last cutpoint")
    }
    stop(
      "'",
      name,
      "' must be a finite value greater than the last value in '",
      cutpoints_name,
      "'"
    )
  }

  invisible(TRUE)
}

#' @title Validate an administrative censoring time
#'
#' @description Checks that an optional administrative censoring time is a
#'   single, finite, strictly positive numeric value after the final cutpoint.
#'
#' @param maxtime `NULL`, or a numeric value giving the administrative censoring
#'   time.
#' @param cutpoints `NULL`, or a numeric vector of interior cutpoints.
#'
#' @noRd
validate_maxtime <- function(maxtime, cutpoints) {
  if (is.null(maxtime)) {
    return(invisible(TRUE))
  }

  validate_endpoint_time(maxtime, cutpoints, "maxtime")
}

#' @title Validate a piecewise hazard vector
#'
#' @description Checks that a finite non-negative hazard vector has one value
#'   for every piecewise interval.
#'
#' @param hazard A numeric vector of piecewise-constant hazard rates.
#' @param cutpoints `NULL`, or a numeric vector of interior cutpoints.
#' @param name A single character string naming the hazard argument in error
#'   messages. The default is `"hazard"`.
#' @param cutpoints_name A single character string naming the cutpoint argument
#'   in error messages. The default is `"cutpoints"`.
#'
#' @noRd
validate_piecewise_hazard <- function(
  hazard,
  cutpoints,
  name = "hazard",
  cutpoints_name = "cutpoints"
) {
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

  if (length(hazard) != length(cutpoints) + 1L) {
    stop(
      "Length of '",
      name,
      "' must be one greater than length of '",
      cutpoints_name,
      "'"
    )
  }

  invisible(TRUE)
}

#' @title Validate a matrix of piecewise hazards
#'
#' @description Checks that posterior hazard draws form a non-empty finite
#'   non-negative matrix with one column per piecewise interval.
#'
#' @param hazard A numeric matrix of posterior piecewise-hazard draws.
#' @param cutpoints `NULL`, or a numeric vector of interior cutpoints.
#' @param name A single character string naming the hazard argument in error
#'   messages. The default is `"hazard"`.
#'
#' @noRd
validate_hazard_matrix <- function(hazard, cutpoints, name = "hazard") {
  if (
    !is.matrix(hazard) ||
      !is.numeric(hazard) ||
      nrow(hazard) == 0 ||
      any(is.na(hazard)) ||
      any(!is.finite(hazard)) ||
      any(hazard < 0)
  ) {
    stop(
      "'",
      name,
      "' must be a non-empty matrix of finite non-negative hazard rates"
    )
  }

  if (ncol(hazard) != length(cutpoints) + 1L) {
    stop(
      "The number of hazard columns must be one greater than length of 'cutpoints'"
    )
  }

  invisible(TRUE)
}

#' @title Validate an enrollment schedule
#'
#' @description Checks the target sample size, piecewise enrollment rates, and
#'   their chronologically ordered knots before simulating accrual.
#'
#' @param lambda A numeric vector of enrollment rates.
#' @param lambda_time `NULL`, or a numeric vector of enrollment-rate change
#'   times.
#' @param N_total A positive integer giving the target sample size.
#'
#' @noRd
validate_enrollment_schedule <- function(lambda, lambda_time, N_total) {
  validate_positive_integer_scalar(N_total, "N_total")

  # A non-positive rate does not define an invertible cumulative intensity and
  # could prevent the process from reaching its target sample size.
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

  if (length(lambda) != length(lambda_time) + 1L) {
    stop(
      "Length of 'lambda' must be one greater than length of 'lambda_time'"
    )
  }

  if (is.null(lambda_time)) {
    return(invisible(TRUE))
  }

  if (
    !is.numeric(lambda_time) ||
      !is.null(dim(lambda_time)) ||
      length(lambda_time) == 0 ||
      any(is.na(lambda_time)) ||
      any(!is.finite(lambda_time)) ||
      any(lambda_time <= 0)
  ) {
    stop(
      "'lambda_time' must be NULL or contain finite positive numeric values"
    )
  }

  if (any(diff(lambda_time) <= 0)) {
    stop("'lambda_time' must be strictly increasing")
  }

  invisible(TRUE)
}

#' @title Validate randomization settings
#'
#' @description Checks that a two-arm block-randomization schedule has valid
#'   block sizes and allocation weights compatible with the target sample size.
#'
#' @param N_total A positive integer giving the target sample size.
#' @param block A positive integer vector of permitted block sizes.
#' @param allocation A length-two positive integer vector giving the control and
#'   treatment allocation weights.
#' @param allocation_name A single character string naming `allocation` in
#'   error messages. The default is `"allocation"`.
#'
#' @noRd
validate_randomization_args <- function(
  N_total,
  block,
  allocation,
  allocation_name = "allocation"
) {
  validate_positive_integer_scalar(N_total, "N_total")

  if (
    !is.numeric(block) ||
      !is.null(dim(block)) ||
      length(block) == 0 ||
      any(is.na(block)) ||
      any(!is.finite(block)) ||
      any(block %% 1 != 0) ||
      any(block <= 0)
  ) {
    stop("'block' must contain positive integer values")
  }

  allocation <- normalize_allocation(allocation, allocation_name)

  if (any(block %% sum(allocation) != 0)) {
    stop(
      "Each 'block' value must be a multiple of sum('",
      allocation_name,
      "') (",
      sum(allocation),
      "); observed 'block' value(s): ",
      paste(block, collapse = ", ")
    )
  }

  if (N_total < sum(block)) {
    stop(
      "'N_total' must be at least sum('block') (",
      sum(block),
      "); observed 'N_total': ",
      N_total
    )
  }

  invisible(allocation)
}

#' @title Validate a null hypothesis value
#'
#' @description Checks that a null value is finite and lies within the support
#'   of probability-scale treatment effects when applicable.
#'
#' @param h0 A single numeric null value.
#' @param method A single character string naming the analysis method.
#' @param single_arm A single logical value indicating whether the design is
#'   single-arm.
#'
#' @noRd
validate_h0 <- function(h0, method, single_arm) {
  if (
    length(h0) != 1 ||
      !is.numeric(h0) ||
      is.na(h0) ||
      !is.finite(h0)
  ) {
    stop("'h0' must be a single finite numeric value")
  }

  if (method == "logrank" && h0 != 0) {
    stop(
      "'h0' must be 0 for log-rank analyses; the supported null is equal ",
      "survival distributions"
    )
  }

  if (
    method %in%
      c("bayes-surv", "bayes-bin", "riskdiff-wald", "riskdiff-fm")
  ) {
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
        if (method %in% c("riskdiff-wald", "riskdiff-fm")) {
          " risk-difference analyses"
        } else {
          " Bayesian analyses"
        }
      )
    }
  }

  invisible(TRUE)
}

#' @title Normalize an interim decision threshold
#'
#' @description Applies the shared scalar-or-one-per-look contract used by
#'   futility, expected-success, and immediate-success thresholds. A scalar is
#'   broadcast to every interim look; any non-scalar vector must have exactly
#'   one value per look.
#'
#' @param threshold `NULL`, or a numeric vector of decision thresholds.
#' @param n_interims A non-negative integer giving the number of interim looks.
#' @param name A single character string naming the threshold in error messages.
#' @param null_disables A single logical value indicating whether `NULL` should
#'   disable the decision rule. The default is `FALSE`.
#'
#' @noRd
normalize_interim_threshold <- function(
  threshold,
  n_interims,
  name,
  null_disables = FALSE
) {
  validate_nonnegative_integer_scalar(n_interims, "n_interims")

  if (n_interims == 0L) {
    return(numeric())
  }

  if (is.null(threshold)) {
    if (null_disables) {
      return(rep.int(0, n_interims))
    }
    stop("'", name, "' must contain finite probabilities in [0, 1]")
  }

  observed <- length(threshold)
  if (!(observed %in% c(1L, n_interims))) {
    direction <- if (observed < n_interims) "too few" else "too many"
    stop(
      "'",
      name,
      "' has ",
      direction,
      " values: expected 1 or ",
      n_interims,
      " (one per interim look), observed ",
      observed
    )
  }

  validate_probability_vector(threshold, name)
  if (observed == 1L) {
    return(rep.int(as.numeric(threshold), n_interims))
  }

  as.numeric(threshold)
}

#' @title Validate ordered interim success thresholds
#'
#' @description Checks that the expected-success boundary does not exceed the
#'   immediate-success boundary at any interim look.
#'
#' @param Sn A normalized numeric vector of expected-success thresholds.
#' @param Qn A normalized numeric vector of immediate-success thresholds.
#'
#' @return `TRUE`, invisibly.
#'
#' @keywords internal
#' @noRd
validate_success_threshold_order <- function(Sn, Qn) {
  if (length(Sn) != length(Qn)) {
    stop("Internal error: 'Sn' and 'Qn' must have the same length")
  }

  invalid <- which(Sn > Qn)
  if (length(invalid) > 0L) {
    stop(
      "'Sn' must be less than or equal to 'Qn' at every interim look; ",
      "Sn > Qn at look(s): ",
      paste(invalid, collapse = ", ")
    )
  }

  invisible(TRUE)
}

#' @title Normalize a deprecated analysis-method name
#'
#' @description Maps the deprecated `"riskdiff"` method to its explicit
#'   `"riskdiff-wald"` replacement and warns once at the public entry point.
#'
#' @param method A single character string naming the analysis method.
#'
#' @return The canonical method name.
#'
#' @noRd
normalize_analysis_method <- function(method) {
  if (identical(method, "riskdiff")) {
    warning(
      "`method = \"riskdiff\"` is deprecated; use ",
      "`method = \"riskdiff-wald\"` instead.",
      call. = FALSE
    )
    return("riskdiff-wald")
  }

  method
}

#' @title Validate analysis-method configuration
#'
#' @description Checks the mutually compatible analysis settings shared by
#'   trial simulation and final analysis.
#'
#' @param method A single character string naming the analysis method.
#' @param alternative A single character string naming the alternative
#'   hypothesis.
#' @param single_arm A single logical value indicating whether the design is
#'   single-arm.
#' @param imputed_final A single logical value indicating whether the final
#'   analysis uses multiply imputed outcomes.
#'
#' @noRd
validate_analysis_configuration <- function(
  method,
  alternative,
  single_arm,
  imputed_final
) {
  if (
    length(alternative) != 1 ||
      !is.character(alternative) ||
      is.na(alternative) ||
      !alternative %in% c("two.sided", "greater", "less")
  ) {
    stop("'alternative' must be one of 'two.sided', 'greater', or 'less'")
  }

  if (
    length(method) != 1 ||
      !is.character(method) ||
      is.na(method) ||
      !method %in%
        c(
          "bayes-surv",
          "bayes-bin",
          "logrank",
          "cox",
          "riskdiff-wald",
          "riskdiff-fm"
        )
  ) {
    stop(
      "'method' must be one of 'bayes-surv', 'bayes-bin', 'logrank', 'cox', ",
      "'riskdiff-wald', or 'riskdiff-fm'"
    )
  }

  validate_logical_scalar(single_arm, "single_arm")
  validate_logical_scalar(imputed_final, "imputed_final")

  if (alternative == "two.sided" && method %in% c("bayes-surv", "bayes-bin")) {
    stop(
      "Bayesian tests can only be used with alternative equal to 'greater' or 'less'"
    )
  }

  if (imputed_final && method == "logrank") {
    stop(
      "The 'logrank' method cannot use 'imputed_final = TRUE' because ",
      "there is no supported frequentist pooling rule for multiple ",
      "imputed final datasets"
    )
  }

  if (
    single_arm &&
      method %in% c("logrank", "cox", "riskdiff-wald", "riskdiff-fm")
  ) {
    stop("The selected method can only be used for two-armed trials")
  }

  invisible(TRUE)
}

#' @title Validate interim analysis looks
#'
#' @description Checks that interim analyses occur at strictly increasing,
#'   positive enrollment counts below the target sample size and, for two-arm
#'   block randomization, after at least one complete block.
#'
#' @param interim_look A positive integer vector giving cumulative enrollment at
#'   each interim analysis.
#' @param N_total A positive integer giving the maximum sample size.
#' @param min_look `NULL`, or a positive integer giving the earliest permitted
#'   interim look. The default is `NULL`.
#'
#' @noRd
validate_interim_looks <- function(interim_look, N_total, min_look = NULL) {
  validate_positive_integer_vector(interim_look, "interim_look")

  if (any(interim_look >= N_total)) {
    stop("'interim_look' must contain values strictly less than 'N_total'")
  }

  if (any(diff(interim_look) <= 0)) {
    stop("'interim_look' must be strictly increasing without duplicates")
  }

  if (!is.null(min_look) && any(interim_look < min_look)) {
    stop(
      "Each 'interim_look' must be at least the block size (",
      min_look,
      ") so that both treatment groups are present at every ",
      "interim analysis. Smallest 'interim_look' given: ",
      min(interim_look),
      "."
    )
  }

  invisible(TRUE)
}

#' @title Validate Bayesian binomial analysis settings
#'
#' @description Checks the Beta prior, computational method, and Monte Carlo
#'   sample size before a Bayesian binomial analysis is run.
#'
#' @param prior_bin A length-two numeric vector of positive Beta shape
#'   parameters.
#' @param bin_method A single character string naming the posterior-probability
#'   calculation.
#' @param N_mcmc A positive integer giving the number of Monte Carlo draws.
#'
#' @noRd
validate_bayes_binomial_args <- function(prior_bin, bin_method, N_mcmc) {
  if (
    !is.numeric(prior_bin) ||
      length(prior_bin) != 2 ||
      anyNA(prior_bin) ||
      any(!is.finite(prior_bin)) ||
      any(prior_bin <= 0)
  ) {
    stop("'prior_bin' must contain two positive finite values")
  }
  if (!bin_method %in% c("mc", "normal", "quadrature")) {
    stop("'bin_method' must be one of 'mc', 'normal', or 'quadrature'")
  }
  if (
    length(N_mcmc) != 1 ||
      !is.numeric(N_mcmc) ||
      is.na(N_mcmc) ||
      !is.finite(N_mcmc) ||
      N_mcmc <= 0 ||
      N_mcmc != floor(N_mcmc)
  ) {
    stop("'N_mcmc' must be a single positive integer")
  }

  invisible(TRUE)
}
