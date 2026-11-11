#' @title Bayesian binomial test for complete binary outcomes
#'
#' @description Updates Beta priors for one or two treatment arms and calculates
#'   the posterior probability of the specified binary-outcome hypothesis.
#'
#' @inheritParams analyse_data
#'
#' @return A list with the posterior probability of success (`success`) and
#'   posterior mean treatment effect (`effect`).
#'
#' @noRd
bayes_binomial_test <- function(
  data,
  single_arm,
  alternative,
  h0,
  prior_bin,
  bin_method,
  N_mcmc
) {
  validate_bayes_binomial_args(prior_bin, bin_method, N_mcmc)
  if (alternative == "two.sided") {
    stop(
      "Bayesian binomial analysis can only be used with alternative equal ",
      "to 'greater' or 'less'"
    )
  }

  treatment_event <- data$event[data$treatment == 1]
  if (length(treatment_event) == 0L) {
    stop("Bayesian binomial analysis requires at least one subject per arm")
  }
  control_event <- if (single_arm) {
    numeric()
  } else {
    data$event[data$treatment == 0]
  }
  if (!single_arm && length(control_event) == 0L) {
    stop("Bayesian binomial analysis requires at least one subject per arm")
  }

  bayes_binomial_from_counts(
    events_control = sum(control_event),
    n_control = length(control_event),
    events_treatment = sum(treatment_event),
    n_treatment = length(treatment_event),
    single_arm = single_arm,
    alternative = alternative,
    h0 = h0,
    prior_bin = prior_bin,
    bin_method = bin_method,
    N_mcmc = N_mcmc
  )
}

#' Calculate a beta-binomial result from arm totals
#'
#' @description Updates the Beta prior using the number of events and subjects
#'   in each arm, then calculates the posterior probability of the specified
#'   hypothesis and the posterior mean treatment effect. These counts are the
#'   sufficient statistics for the completed binary-endpoint analysis.
#'
#' @param events_control A non-negative integer giving the number of control-arm
#'   events. Ignored in a single-arm analysis.
#' @param n_control A non-negative integer giving the number of control-arm
#'   subjects. Ignored in a single-arm analysis.
#' @param events_treatment A non-negative integer giving the number of
#'   treatment-arm events.
#' @param n_treatment A positive integer giving the number of treatment-arm
#'   subjects.
#' @inheritParams analyse_data
#'
#' @return A list containing the posterior probability of the alternative
#'   (`success`) and posterior mean event probability or treatment-control
#'   difference (`effect`). Monte Carlo analyses also retain the number of
#'   posterior draws satisfying the alternative.
#'
#' @keywords internal
#' @noRd
bayes_binomial_from_counts <- function(
  events_control,
  n_control,
  events_treatment,
  n_treatment,
  single_arm,
  alternative,
  h0,
  prior_bin,
  bin_method,
  N_mcmc
) {
  beta_binomial_stats <- function(events, n) {
    alpha <- prior_bin[1] + events
    beta <- prior_bin[2] + n - events
    list(
      alpha = alpha,
      beta = beta,
      mean = alpha / (alpha + beta),
      variance = (alpha * beta) / ((alpha + beta)^2 * (alpha + beta + 1))
    )
  }

  beta_binomial_difference_success <- function(treatment, control) {
    integrand <- function(x) {
      treatment_density <- dbeta(x, treatment$alpha, treatment$beta)
      control_threshold <- x - h0
      if (alternative == "greater") {
        treatment_density *
          pbeta(control_threshold, control$alpha, control$beta)
      } else {
        treatment_density *
          (1 - pbeta(control_threshold, control$alpha, control$beta))
      }
    }

    # Numerical integration can exceed one by a few machine epsilons in an
    # extreme posterior tail. Keep the result on the probability scale so a
    # threshold of one remains an impossible strict crossing.
    pmin(1, pmax(0, integrate(integrand, lower = 0, upper = 1)$value))
  }

  treatment_stats <- beta_binomial_stats(events_treatment, n_treatment)

  if (single_arm) {
    if (bin_method == "mc") {
      effect_draws <- rbeta(
        N_mcmc,
        treatment_stats$alpha,
        treatment_stats$beta
      )
      success_indicator <- if (alternative == "greater") {
        effect_draws > h0
      } else {
        effect_draws < h0
      }
      success <- mean(success_indicator)
      effect <- mean(effect_draws)
    } else {
      effect <- treatment_stats$mean
      if (bin_method == "normal") {
        effect_se <- sqrt(treatment_stats$variance)
        success <- if (alternative == "greater") {
          1 - pnorm(h0, mean = effect, sd = effect_se)
        } else {
          pnorm(h0, mean = effect, sd = effect_se)
        }
      } else {
        success <- if (alternative == "greater") {
          1 - pbeta(h0, treatment_stats$alpha, treatment_stats$beta)
        } else {
          pbeta(h0, treatment_stats$alpha, treatment_stats$beta)
        }
      }
    }
    result <- list(success = success, effect = effect)
    if (bin_method == "mc") {
      result <- set_analysis_mc_counts(
        result,
        successes = sum(success_indicator),
        draws = length(success_indicator)
      )
    }
    return(result)
  }

  control_stats <- beta_binomial_stats(events_control, n_control)

  if (bin_method == "mc") {
    effect_draws <-
      rbeta(N_mcmc, treatment_stats$alpha, treatment_stats$beta) -
      rbeta(N_mcmc, control_stats$alpha, control_stats$beta)
    success_indicator <- if (alternative == "greater") {
      effect_draws > h0
    } else {
      effect_draws < h0
    }
    success <- mean(success_indicator)
    effect <- mean(effect_draws)
  } else if (bin_method == "normal") {
    effect <- treatment_stats$mean - control_stats$mean
    effect_se <- sqrt(treatment_stats$variance + control_stats$variance)
    if (alternative == "greater") {
      success <- 1 - pnorm(h0, mean = effect, sd = effect_se)
    } else if (alternative == "less") {
      success <- pnorm(h0, mean = effect, sd = effect_se)
    }
  } else if (bin_method == "quadrature") {
    effect <- treatment_stats$mean - control_stats$mean
    success <- beta_binomial_difference_success(treatment_stats, control_stats)
    check_beta_quadrature(success, treatment_stats, control_stats, h0, alternative)
  }

  result <- list(success = success, effect = effect)
  if (bin_method == "mc") {
    result <- set_analysis_mc_counts(
      result,
      successes = sum(success_indicator),
      draws = length(success_indicator)
    )
  }
  result
}

# Numerical checks for two-arm Beta quadrature. These diagnose probabilities;
# they do not replace the integration result or change the stopping rules.
check_beta_quadrature <- function(
  probability,
  treatment,
  control,
  h0,
  alternative
) {
  tolerance <- 1e-6
  symmetry <- h0 == 0 &&
    treatment$alpha == control$alpha &&
    treatment$beta == control$beta &&
    abs(probability - 0.5) > tolerance

  # Cantelli's inequality for D = p_treatment - p_control - h0.
  # With negative mean, Pr(D < 0) >= mean^2 / (variance + mean^2);
  # with positive mean, Pr(D < 0) <= variance / (variance + mean^2).
  mean_difference <- treatment$mean - control$mean - h0
  variance <- treatment$variance + control$variance
  denominator <- variance + mean_difference^2
  moments <- FALSE
  if (
    is.finite(mean_difference) &&
      is.finite(variance) &&
      variance >= 0 &&
      is.finite(denominator) &&
      denominator > 0
  ) {
    probability_less <- if (alternative == "less") {
      probability
    } else {
      1 - probability
    }
    if (mean_difference < 0) {
      moments <- probability_less < mean_difference^2 / denominator - tolerance
    } else if (mean_difference > 0) {
      moments <- probability_less > variance / denominator + tolerance
    }
  }
  if (symmetry || moments) {
    warn_quadrature_diagnostics(c(
      flagged_calculations = 1,
      symmetry = as.numeric(symmetry),
      moments = as.numeric(moments)
    ))
  }
  invisible(NULL)
}

new_quadrature_diagnostics <- function() {
  diagnostics <- new.env(parent = emptyenv())
  diagnostics$counts <- c(flagged_calculations = 0, symmetry = 0, moments = 0)
  diagnostics
}

warn_quadrature_diagnostics <- function(counts) {
  if (counts[["flagged_calculations"]] == 0) {
    return(invisible(NULL))
  }
  message <- paste0(
    "Numerical accuracy checks flagged ",
    counts[["flagged_calculations"]],
    " two-arm quadrature probability calculation(s) (symmetry: ",
    counts[["symmetry"]],
    "; posterior mean-variance bounds: ",
    counts[["moments"]],
    "). Success and stopping decisions may be unreliable. ",
    "Check these probabilities using sufficiently precise Monte Carlo Beta ",
    "draws (bin_method = \"mc\"). The quadrature results have not been changed."
  )
  warning(structure(
    list(message = message, call = NULL, counts = counts),
    class = c("goldilocks_quadrature_warning", "warning", "condition")
  ))
}

# Catch only this diagnostic warning; leave other warning policies unchanged.
with_quadrature_warnings <- function(expr, diagnostics) {
  withCallingHandlers(expr, goldilocks_quadrature_warning = function(w) {
    diagnostics$counts <- diagnostics$counts + w$counts
    invokeRestart("muffleWarning")
  })
}

finish_quadrature_warnings <- function(result, diagnostics, warn = TRUE) {
  if (diagnostics$counts[["flagged_calculations"]] > 0) {
    attr(result, "quadrature_diagnostics") <- c(
      diagnostics$counts,
      absolute_tolerance = 1e-6
    )
    if (warn) warn_quadrature_diagnostics(diagnostics$counts)
  }
  result
}

# Workers return counts to the calling process, which issues the run warning.
capture_quadrature_warnings <- function(expr) {
  diagnostics <- new_quadrature_diagnostics()
  value <- with_quadrature_warnings(expr, diagnostics)
  finish_quadrature_warnings(value, diagnostics, warn = FALSE)
}

combine_quadrature_warnings <- function(results) {
  diagnostics <- new_quadrature_diagnostics()
  for (result in results) {
    counts <- attr(result, "quadrature_diagnostics", exact = TRUE)
    if (!is.null(counts)) {
      diagnostics$counts <- diagnostics$counts +
        counts[names(diagnostics$counts)]
    }
  }
  diagnostics
}
