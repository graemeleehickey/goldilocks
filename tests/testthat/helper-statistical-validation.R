mc_binomial_acceptance_interval <- function(
  target,
  n,
  level = 0.9999
) {
  tail_probability <- (1 - level) / 2
  stats::qbinom(
    c(tail_probability, 1 - tail_probability),
    size = n,
    prob = target
  ) /
    n
}

expect_mc_rate <- function(
  indicator,
  target,
  estimand,
  level = 0.9999
) {
  n <- length(indicator)
  estimate <- mean(indicator)
  mcse <- sqrt(target * (1 - target) / n)
  limits <- mc_binomial_acceptance_interval(target, n, level)

  testthat::expect_true(
    estimate >= limits[[1L]] && estimate <= limits[[2L]],
    info = sprintf(
      paste0(
        "%s: target = %.6f; estimate = %.6f; MCSE = %.6f; ",
        "acceptance interval = [%.6f, %.6f]"
      ),
      estimand,
      target,
      estimate,
      mcse,
      limits[[1L]],
      limits[[2L]]
    )
  )

  invisible(list(
    estimand = estimand,
    target = target,
    estimate = estimate,
    mcse = mcse,
    lower = limits[[1L]],
    upper = limits[[2L]],
    n = n
  ))
}

expect_mc_close <- function(
  estimate,
  target,
  mcse,
  estimand,
  reference_mcse = 0,
  standard_errors = 4
) {
  combined_mcse <- sqrt(mcse^2 + reference_mcse^2)
  tolerance <- standard_errors * combined_mcse

  testthat::expect_true(
    abs(estimate - target) <= tolerance,
    info = sprintf(
      paste0(
        "%s: target = %.6f; estimate = %.6f; MCSE = %.6f; ",
        "tolerance = %.6f"
      ),
      estimand,
      target,
      estimate,
      combined_mcse,
      tolerance
    )
  )

  invisible(list(
    estimand = estimand,
    target = target,
    estimate = estimate,
    mcse = combined_mcse,
    tolerance = tolerance
  ))
}

gamma_variance_mcse <- function(shape, rate, n) {
  variance <- shape / rate^2
  fourth_central_moment <- (3 + 6 / shape) * variance^2

  sqrt(
    (fourth_central_moment -
      ((n - 3) / (n - 1)) * variance^2) /
      n
  )
}
