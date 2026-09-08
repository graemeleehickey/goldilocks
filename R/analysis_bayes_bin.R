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

    integrate(integrand, lower = 0, upper = 1)$value
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
