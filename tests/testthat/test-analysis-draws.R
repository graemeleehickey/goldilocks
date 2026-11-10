draw_count_cases <- function(method, bin_method) {
  common <- list(
    hazard_treatment = 0.05, hazard_control = 0.08,
    N_total = 24, lambda = 4, alternative = "less",
    method = method, bin_method = bin_method, N_impute = 3, N_mcmc = 20
  )
  fixed <- c(common, list(interim_look = 12, end_of_study = 12))
  horizon <- c(common, list(
    endpoint_schedule = data.frame(N = c(12, 18), horizon = c(8, 12)),
    final_horizon = 12
  ))
  repeated <- list(N_trials = 1, seed = 8301, backend = "sequential")
  rules <- list(Fn = 0, Sn = 1, Qn = 1, return_trace = TRUE)
  trial_results <- function(x) list(summary = lapply(x$summary, identity), trace = x$trace)
  simulation_results <- function(x) list(sims = lapply(x$sims, identity), traces = x$traces)
  path_results <- function(x) list(trials = x$trials, failures = x$failures)

  list(
    survival_adapt = list(fun = "survival_adapt", args = c(fixed, rules), result = trial_results),
    sim_trials = list(fun = "sim_trials", args = c(fixed, repeated, rules), result = simulation_results),
    sim_trial_paths = list(fun = "sim_trial_paths", args = c(fixed, repeated), result = path_results),
    sim_horizon_trials = list(fun = "sim_horizon_trials", args = c(horizon, repeated, rules), result = simulation_results),
    sim_horizon_paths = list(fun = "sim_horizon_paths", args = c(horizon, repeated), result = path_results),
    evaluate_interim = list(
      fun = "evaluate_interim",
      args = list(
        data = data.frame(
          id = 1:8, treatment = rep(0:1, 4), enrollment = 0:7,
          time = 8:1, event = c(1, 0, 0, 1, 0, 0, 0, 0),
          status = c("event", "pending", "censored", "event", rep("pending", 4))
        ),
        data_cut = 8, look = 1, N_total = 12, end_of_study = 10,
        alternative = "less", method = method, bin_method = bin_method,
        N_impute = 3, N_mcmc = 20, seed = 8301
      ),
      result = function(x) x[c("decision", "probabilities", "monte_carlo", "diagnostics", "trace")]
    )
  )
}

test_that("deterministic analyses ignore N_mcmc across trial functions", {
  for (bin_method in c("normal", "quadrature")) {
    cases <- draw_count_cases("bayes-bin", bin_method)
    for (name in names(cases)) {
      case <- cases[[name]]
      set.seed(8301)
      reference <- case$result(do.call(case$fun, case$args))
      for (draws in list(NULL, 0, NA_real_, 1.5, "unused")) {
        args <- case$args
        args["N_mcmc"] <- list(draws)
        set.seed(8301)
        result <- do.call(case$fun, args)
        expect_identical(case$result(result), reference, info = paste(name, bin_method))
      }
    }
  }
})

test_that("frequentist analyses accept NULL for an unused draw count", {
  for (case in draw_count_cases("logrank", "mc")) {
    set.seed(8301)
    reference <- case$result(do.call(case$fun, case$args))
    case$args["N_mcmc"] <- list(NULL)
    set.seed(8301)
    expect_identical(case$result(do.call(case$fun, case$args)), reference)
  }
})

test_that("posterior sampling still requires a positive integer draw count", {
  cases <- c(
    draw_count_cases("bayes-bin", "mc"),
    draw_count_cases("bayes-surv", "quadrature")
  )
  for (case in cases) {
    for (draws in list(NULL, 0, 1.5, NA_real_)) {
      args <- case$args
      args["N_mcmc"] <- list(draws)
      expect_error(do.call(case$fun, args), "N_mcmc.*positive integer")
    }
  }
})

test_that("deterministic binomial calculations do not require N_mcmc", {
  for (single_arm in c(FALSE, TRUE)) {
    data <- data.frame(
      treatment = rep(0:1, each = 6),
      event = c(1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0)
    )
    if (single_arm) data <- data[data$treatment == 1, ]
    for (bin_method in c("normal", "quadrature")) {
      args <- list(
        data = data, single_arm = single_arm, alternative = "less",
        h0 = if (single_arm) 0.5 else 0, prior_bin = c(1, 1), bin_method = bin_method
      )
      without_draws <- do.call(bayes_binomial_test, args)
      expect_identical(
        without_draws,
        do.call(bayes_binomial_test, c(args, list(N_mcmc = 100)))
      )
    }
  }
})

test_that("ignored draw counts remain available in interim settings", {
  case <- draw_count_cases("bayes-bin", "quadrature")$evaluate_interim
  for (draws in list(NULL, 0, "unused")) {
    case$args["N_mcmc"] <- list(draws)
    result <- do.call(case$fun, case$args)
    expect_true("N_mcmc" %in% names(result$metadata$design))
    expect_identical(result$metadata$design$N_mcmc, draws)
  }
})
