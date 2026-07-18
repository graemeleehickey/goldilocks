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

#' @title Validate a positive integer vector
#'
#' @description Checks that an input vector contains only finite, strictly
#'   positive integers.
#'
#' @noRd
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

#' @title Validate a Gamma prior
#'
#' @description Checks that a Gamma prior supplies two finite, strictly positive
#'   shape and rate parameters.
#'
#' @noRd
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

#' @title Validate piecewise cutpoints
#'
#' @description Checks that a piecewise model has finite, strictly increasing
#'   cutpoints beginning at zero.
#'
#' @noRd
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

#' @title Validate an endpoint time
#'
#' @description Checks that an analysis endpoint is finite and positive and,
#'   when required, lies after the final piecewise cutpoint.
#'
#' @noRd
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

#' @title Validate an administrative censoring time
#'
#' @description Checks that an optional administrative censoring time is a
#'   single, finite, strictly positive numeric value.
#'
#' @noRd
validate_maxtime <- function(maxtime) {
  if (
    !is.null(maxtime) &&
      (
        length(maxtime) != 1 ||
          !is.numeric(maxtime) ||
          is.na(maxtime) ||
          !is.finite(maxtime) ||
          maxtime <= 0
      )
  ) {
    stop("'maxtime' must be a single finite positive value or NULL")
  }

  invisible(TRUE)
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

  if (length(hazard) != length(cutpoints)) {
    stop("Length of 'cutpoints' must be equal to length of '", name, "'")
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
    stop("'", name, "' must be a non-empty matrix of finite non-negative hazard rates")
  }

  if (ncol(hazard) != length(cutpoints)) {
    stop("The length of the hazard rates and cutpoints do not match")
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
    stop("Sum of 'allocation' must be a multiple of 'block'")
  }

  if (N_total < sum(block)) {
    stop("Number of subjects must be at least the size of 'block'")
  }

  invisible(TRUE)
}

#' @title Validate a null hypothesis value
#'
#' @description Checks that a null value is finite and lies within the support
#'   of Bayesian probability-scale treatment effects when applicable.
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
      !method %in% c("bayes-surv", "bayes-bin", "logrank", "cox", "chisq")
  ) {
    stop("'method' must be one of 'bayes-surv', 'bayes-bin', 'logrank', 'cox', or 'chisq'")
  }

  validate_logical_scalar(single_arm, "single_arm")
  validate_logical_scalar(imputed_final, "imputed_final")

  if (alternative == "two.sided" && method %in% c("bayes-surv", "bayes-bin")) {
    stop(
      "Bayesian tests can only be used with alternative equal to 'greater' or 'less'"
    )
  }

  if (alternative != "two.sided" && method == "chisq") {
    stop("The chi-square test can only be applied as a two-sided test")
  }

  if (imputed_final && method %in% c("logrank", "chisq")) {
    stop(
      "The 'logrank' and 'chisq' methods cannot use ",
      "'imputed_final = TRUE' because ",
      "there is no supported frequentist pooling rule for multiple ",
      "imputed final datasets"
    )
  }

  if (single_arm && method %in% c("logrank", "cox", "chisq")) {
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
validate_bayes_binomial_args <- function(bin_prior, bin_method, N_mcmc) {
  if (
    length(bin_prior) != 2 ||
      any(!is.finite(bin_prior)) ||
      any(bin_prior <= 0)
  ) {
    stop("'bin_prior' must contain two positive finite values")
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
