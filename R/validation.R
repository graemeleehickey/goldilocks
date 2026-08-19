#' @title Validate one probability
#'
#' @description Checks that an input is a single finite probability, optionally
#'   excluding one to support distributions with an open upper bound.
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

#' Normalize arm-specific loss-to-follow-up proportions
#'
#' @description Converts a shared scalar or a named arm-specific vector into
#'   the canonical treatment-group order used by the simulation. Requiring
#'   names for two-arm vectors prevents positional confusion between APIs that
#'   list treatment and control in different orders.
#'
#' @param prop_loss One shared probability or a named arm-specific vector.
#' @param single_arm Whether the simulation has only a treatment arm.
#'
#' @return A named probability vector containing `treatment` for a single-arm
#'   design or `control` followed by `treatment` for a two-arm design.
#'
#' @keywords internal
#' @noRd
normalize_prop_loss <- function(prop_loss, single_arm) {
  validate_probability_vector(prop_loss, "prop_loss")

  if (single_arm) {
    if (length(prop_loss) != 1L) {
      stop("'prop_loss' must be a single probability for a single-arm design")
    }
    return(c(treatment = unname(prop_loss)))
  }

  arm_names <- c("control", "treatment")
  if (length(prop_loss) == 1L) {
    normalized <- rep(unname(prop_loss), 2L)
    names(normalized) <- arm_names
    return(normalized)
  }

  supplied_names <- names(prop_loss)
  valid_arm_vector <- length(prop_loss) == 2L &&
    !is.null(supplied_names) &&
    !anyNA(supplied_names) &&
    !anyDuplicated(supplied_names) &&
    setequal(supplied_names, arm_names)
  if (!valid_arm_vector) {
    stop(
      "'prop_loss' must be a single probability or a length-two vector ",
      "named 'control' and 'treatment' for a two-arm design"
    )
  }

  normalized <- unname(prop_loss[arm_names])
  names(normalized) <- arm_names
  normalized
}

#' @title Validate one positive integer
#'
#' @description Checks that an input is a single finite, strictly positive
#'   integer suitable for a sample size or count.
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
#' @description Checks that a Gamma prior supplies finite, strictly positive
#'   shape and rate parameters, then broadcasts a length-two vector over all
#'   piecewise intervals.
#'
#' @noRd
normalize_gamma_prior <- function(
  prior_surv,
  n_intervals,
  name = "prior_surv"
) {
  valid_values <- is.numeric(prior_surv) &&
    !anyNA(prior_surv) &&
    all(is.finite(prior_surv)) &&
    all(prior_surv > 0)
  is_vector <- is.null(dim(prior_surv)) && length(prior_surv) == 2L
  is_interval_matrix <- is.matrix(prior_surv) &&
    identical(dim(prior_surv), c(2L, as.integer(n_intervals)))

  if (!valid_values || (!is_vector && !is_interval_matrix)) {
    stop(
      "'",
      name,
      "' must be a length-two positive finite numeric vector or a 2 x ",
      n_intervals,
      " matrix with shape in row 1 and rate in row 2"
    )
  }

  if (is_vector) {
    prior_surv <- matrix(
      rep(prior_surv, n_intervals),
      nrow = 2L,
      dimnames = list(c("shape", "rate"), NULL)
    )
  }

  prior_surv
}

#' @title Validate piecewise cutpoints
#'
#' @description Checks that a piecewise model has finite, positive, strictly
#'   increasing interior cutpoints, or no cutpoints for a constant hazard.
#'
#' @noRd
validate_cutpoints <- function(cutpoints) {
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
    stop("'cutpoints' must be NULL or contain finite positive numeric values")
  }

  # A non-increasing knot vector creates zero-width or backward intervals,
  # which makes piecewise exposure and event allocation undefined.
  if (any(diff(cutpoints) <= 0)) {
    stop("'cutpoints' must be strictly increasing")
  }

  invisible(TRUE)
}

#' @title Validate an endpoint time
#'
#' @description Checks that an analysis endpoint is finite and positive and
#'   lies after the final piecewise cutpoint.
#'
#' @noRd
validate_endpoint_time <- function(endpoint, cutpoints, name) {
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
    stop("'", name, "' must be a finite value greater than the last cutpoint")
  }

  invisible(TRUE)
}

#' @title Validate an administrative censoring time
#'
#' @description Checks that an optional administrative censoring time is a
#'   single, finite, strictly positive numeric value after the final cutpoint.
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
#' @noRd
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

  if (length(hazard) != length(cutpoints) + 1L) {
    stop(
      "Length of '",
      name,
      "' must be one greater than length of 'cutpoints'"
    )
  }

  invisible(TRUE)
}

#' @title Validate a matrix of piecewise hazards
#'
#' @description Checks that posterior hazard draws form a non-empty finite
#'   non-negative matrix with one column per piecewise interval.
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
#' @noRd
validate_randomization_args <- function(N_total, block, allocation) {
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

  if (length(allocation) != 2) {
    stop("'allocation' must contain two positive integer values")
  }

  if (
    !is.numeric(allocation) ||
      !is.null(dim(allocation)) ||
      any(is.na(allocation)) ||
      any(!is.finite(allocation)) ||
      any(allocation %% 1 != 0)
  ) {
    stop("All values of 'allocation' must be integer values")
  }

  if (any(allocation <= 0)) {
    stop("'allocation' must contain two positive integer values")
  }

  if (any(block %% sum(allocation) != 0)) {
    stop(
      "Each 'block' value must be a multiple of sum('allocation') (",
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

  invisible(TRUE)
}

#' @title Validate a null hypothesis value
#'
#' @description Checks that a null value is finite and lies within the support
#'   of probability-scale treatment effects when applicable.
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

  if (method %in% c("bayes-surv", "bayes-bin", "riskdiff")) {
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
        if (method == "riskdiff") {
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
#'   futility and expected-success thresholds. A scalar is broadcast to every
#'   interim look; any non-scalar vector must have exactly one value per look.
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

#' @title Validate analysis-method configuration
#'
#' @description Checks the mutually compatible analysis settings shared by
#'   trial simulation and final analysis.
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
          "riskdiff"
        )
  ) {
    stop(
      "'method' must be one of 'bayes-surv', 'bayes-bin', 'logrank', 'cox', ",
      "or 'riskdiff'"
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

  if (single_arm && method %in% c("logrank", "cox", "riskdiff")) {
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
