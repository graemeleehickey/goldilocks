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
