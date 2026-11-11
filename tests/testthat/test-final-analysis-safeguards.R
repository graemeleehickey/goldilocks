final_safeguard_data <- function() {
  data.frame(
    treatment = rep(0:1, each = 10),
    time = c(0.1, 0.3, 0.5, rep(1, 7), 0.2, 0.4, 0.6, 0.8, rep(1, 6)),
    event = c(rep(1, 3), rep(0, 7), rep(1, 4), rep(0, 6)),
    loss_to_fu = FALSE,
    subject_impute_success = FALSE
  )
}

final_safeguard_args <- function(data = final_safeguard_data()) {
  list(
    data_in = data,
    cutpoints = c(0.4, 0.7),
    prior_surv_final = c(0.1, 0.1),
    N_mcmc = 100,
    single_arm = FALSE,
    imputed_final = FALSE,
    method = "riskdiff-fm",
    N_impute = 1,
    alternative = "less",
    h0 = 0.1,
    prior_bin = c(1, 1),
    bin_method = "mc",
    binary_imputation = "event-time",
    empty_interval = "prior",
    end_of_study = 1,
    rmst_tau = 0.8
  )
}

test_that("complete final outcomes use the same selected test with either flag", {
  for (method in c(
    "riskdiff-fm",
    "riskdiff-wald",
    "cox",
    "rmst",
    "bayes-surv",
    "bayes-bin"
  )) {
    args <- final_safeguard_args()
    args$method <- method
    set.seed(9271)
    direct <- do.call(analyse_final, args)
    direct_rng <- .Random.seed
    args$imputed_final <- TRUE
    set.seed(9271)
    flagged <- do.call(analyse_final, args)
    expect_identical(flagged, direct, info = method)
    expect_identical(.Random.seed, direct_rng, info = method)
  }
})

test_that("conjugate survival finals ignore the imputation flag and count", {
  local_mocked_bindings(
    impute_data = function(...) stop("Final imputation should not be reached")
  )
  for (single_arm in c(FALSE, TRUE)) {
    for (empty in c(FALSE, TRUE)) {
      data <- final_safeguard_data()
      data$time[c(4, 15)] <- c(0.25, 0.55)
      if (empty) {
        data$event[data$time > 0.65] <- 0
        data$time <- pmin(data$time, 0.65)
      }
      data$loss_to_fu <- data$event == 0 & data$time < 1
      data$subject_impute_success <- data$loss_to_fu
      if (single_arm) {
        data <- data[data$treatment == 1, ]
      }
      args <- final_safeguard_args(data)
      args$method <- "bayes-surv"
      args$single_arm <- single_arm
      args$h0 <- if (single_arm) 0.5 else -0.1
      args$prior_surv_final <- list(
        treatment = matrix(c(2, 3, 3, 4, 4, 5), nrow = 2)
      )
      if (!single_arm) {
        args$prior_surv_final$control <- matrix(c(1, 2, 2, 3, 3, 4), nrow = 2)
      }
      for (alternative in c("less", "greater")) {
        args$alternative <- alternative
        args$imputed_final <- FALSE
        set.seed(9274)
        direct <- do.call(analyse_final, args)
        direct_rng <- .Random.seed
        args$imputed_final <- TRUE
        for (n_impute in c(1, 7)) {
          args$N_impute <- n_impute
          set.seed(9274)
          expect_identical(do.call(analyse_final, args), direct)
          expect_identical(.Random.seed, direct_rng)
        }
      }
    }
  }
})

test_that("direct survival finals agree with analytical censored-data inference", {
  data <- data.frame(
    treatment = 1,
    time = c(0.2, 0.4, 0.6, 1),
    event = c(1, 1, 0, 0),
    loss_to_fu = c(FALSE, FALSE, TRUE, FALSE),
    subject_impute_success = c(FALSE, FALSE, TRUE, FALSE)
  )
  args <- final_safeguard_args(data)
  args["cutpoints"] <- list(NULL)
  args$method <- "bayes-surv"
  args$single_arm <- TRUE
  args$imputed_final <- TRUE
  args$prior_surv_final <- c(1, 2)
  args$N_mcmc <- 100000
  args$h0 <- 0.5

  # Two observed events and 2.2 units of exposure, including the censored
  # participant's 0.6 units, give a Gamma(3, 4.2) posterior.
  shape <- 3
  rate <- 4.2
  limit <- -log1p(-args$h0)
  probability <- pgamma(limit, shape, rate = rate)
  mean_risk <- 1 - (rate / (rate + 1))^shape
  variance_risk <- (rate / (rate + 2))^shape - (1 - mean_risk)^2
  set.seed(9275)
  result <- do.call(analyse_final, args)
  expect_lt(
    abs(result[1] - probability),
    6 * sqrt(probability * (1 - probability) / args$N_mcmc)
  )
  expect_lt(abs(result[2] - mean_risk), 6 * sqrt(variance_risk / args$N_mcmc))

  # Integrate over the one missing outcome: an event during the remaining
  # 0.4 units, or survival to the endpoint. Both completed-data summaries
  # average back to the observed-data summaries without simulation.
  remaining <- 0.4
  density <- function(u) shape * rate^shape / (rate + u)^(shape + 1)
  survival <- (rate / (rate + remaining))^shape
  averaged_probability <- integrate(
    function(u) {
      pgamma(limit, shape + 1, rate = rate + u) * density(u)
    },
    0,
    remaining
  )$value +
    survival * pgamma(limit, shape, rate = rate + remaining)
  averaged_risk <- integrate(
    function(u) {
      (1 - ((rate + u) / (rate + u + 1))^(shape + 1)) * density(u)
    },
    0,
    remaining
  )$value +
    survival * (1 - ((rate + remaining) / (rate + remaining + 1))^shape)
  expect_equal(averaged_probability, probability, tolerance = 1e-10)
  expect_equal(averaged_risk, mean_risk, tolerance = 1e-10)
})

test_that("other survival interval policies and binary finals still impute", {
  data <- final_safeguard_data()
  data$time[4] <- 0.25
  data$loss_to_fu[4] <- TRUE
  data$subject_impute_success[4] <- TRUE
  args <- final_safeguard_args(data)
  args$imputed_final <- TRUE
  args$N_impute <- 3
  args$prior_surv_final <- normalize_gamma_prior(
    args$prior_surv_final,
    3L,
    FALSE,
    "prior_surv_final"
  )
  calls <- 0L
  original_impute <- impute_data
  local_mocked_bindings(impute_data = function(...) {
    calls <<- calls + 1L
    original_impute(...)
  })
  for (policy in c("propagate", "error")) {
    args$method <- "bayes-surv"
    args$empty_interval <- policy
    calls <- 0L
    set.seed(9276)
    result <- do.call(analyse_final, args)
    expect_true(all(is.finite(result)))
    expect_identical(calls, 3L)
  }
  args$method <- "bayes-bin"
  args$empty_interval <- "prior"
  calls <- 0L
  set.seed(9276)
  expect_true(all(is.finite(do.call(analyse_final, args))))
  expect_identical(calls, 3L)

  args$method <- "bayes-surv"
  args$empty_interval <- "error"
  args$data_in$event[args$data_in$time > 0.65] <- 0
  args$data_in$time <- pmin(args$data_in$time, 0.65)
  calls <- 0L
  expect_error(do.call(analyse_final, args), "zero subjects")
  expect_identical(calls, 0L)
})

test_that("survival final flags preserve trial predictions and final results", {
  args <- list(
    hazard_treatment = c(0.03, 0.05),
    hazard_control = c(0.04, 0.06),
    cutpoints = 6,
    N_total = 60,
    lambda = 5,
    interim_look = c(30, 45),
    end_of_study = 12,
    prop_loss = 0.35,
    prior_surv = c(0.1, 0.1),
    prior_surv_final = c(1, 5),
    method = "bayes-surv",
    alternative = "less",
    N_impute = 4,
    N_mcmc = 50,
    Fn = 0.01,
    Sn = 0.99,
    return_trace = TRUE
  )
  set.seed(9277)
  direct <- do.call(survival_adapt, c(args, list(imputed_final = FALSE)))
  direct_rng <- .Random.seed
  set.seed(9277)
  flagged <- do.call(survival_adapt, c(args, list(imputed_final = TRUE)))
  expect_gt(nrow(direct$trace), 0L)
  expect_true(is.finite(direct$summary$post_prob_ha))
  expect_identical(flagged$trace, direct$trace)
  expect_identical(flagged$posterior_diagnostics, direct$posterior_diagnostics)
  columns <- setdiff(names(direct$summary), "imputed_final")
  expect_identical(flagged$summary[columns], direct$summary[columns])
  expect_identical(.Random.seed, direct_rng)
})

test_that("complete FM boundary tables retain the score test and its null variance", {
  args <- final_safeguard_args()
  args$data_in$time <- 1
  args$data_in$event <- 0
  args$imputed_final <- TRUE
  result <- do.call(analyse_final, args)
  # At zero events and a positive margin, the null-constrained variance is
  # h0 * (1 - h0) / n_treatment, not the zero plug-in variance.
  expected <- pnorm(sqrt(10 * 0.1 / (1 - 0.1)))
  expect_equal(result, c(expected, 0), tolerance = 1e-7)
  expect_lt(result[1], 0.975)

  for (event in list(rep(0, 20), rep(1, 20), final_safeguard_data()$event)) {
    args$data_in$event <- event
    for (alternative in c("less", "greater", "two.sided")) {
      args$alternative <- alternative
      for (h0 in c(-0.1, 0, 0.1)) {
        args$h0 <- h0
        args$imputed_final <- FALSE
        direct <- do.call(analyse_final, args)
        args$imputed_final <- TRUE
        expect_identical(do.call(analyse_final, args), direct)
      }
    }
  }
})

test_that("genuinely missing final outcomes cannot silently switch FM to Wald", {
  data <- final_safeguard_data()
  data$time[4] <- 0.25
  data$loss_to_fu[4] <- TRUE
  data$subject_impute_success[4] <- TRUE
  args <- final_safeguard_args(data)
  args$imputed_final <- TRUE
  args$N_impute <- 5
  set.seed(9272)
  rng <- .Random.seed
  expect_error(do.call(analyse_final, args), "riskdiff-fm.*pooling rule")
  expect_identical(.Random.seed, rng)

  # Observed complete-case FM analysis remains supported when requested.
  args$imputed_final <- FALSE
  actual <- do.call(analyse_final, args)
  expected <- risk_difference_fm_test_checked(
    data[!data$loss_to_fu, ],
    1,
    "less",
    0.1
  )
  expect_equal(actual, c(expected$success, expected$estimate))
})

test_that("FM imputation with dropout is rejected before simulation setup", {
  local_mocked_bindings(
    sim_comp_data = function(...) stop("data generation was reached"),
    resolve_sim_execution = function(...) stop("worker setup was reached")
  )
  args <- list(
    hazard_control = 0.2,
    hazard_treatment = 0.1,
    N_total = 40,
    end_of_study = 1,
    method = "riskdiff-fm",
    imputed_final = TRUE,
    N_impute = 5
  )
  for (dropout in list(
    0.1,
    c(control = 0, treatment = 0.1),
    c(control = 0.1, treatment = 0)
  )) {
    args$prop_loss <- dropout
    expect_error(do.call(survival_adapt, args), "riskdiff-fm.*pooling rule")
    expect_error(do.call(sim_trials, args), "riskdiff-fm.*pooling rule")
  }
})

test_that("FM simulations without dropout preserve results and output columns", {
  args <- list(
    hazard_control = 0.2,
    hazard_treatment = 0.1,
    N_total = 40,
    end_of_study = 1,
    prop_loss = 0,
    method = "riskdiff-fm",
    alternative = "less",
    h0 = 0.1,
    N_impute = 1,
    N_trials = 3,
    backend = "sequential",
    seed = 9273
  )
  direct <- do.call(sim_trials, c(args, list(imputed_final = FALSE)))
  flagged <- do.call(sim_trials, c(args, list(imputed_final = TRUE)))
  expect_identical(names(flagged), names(direct))
  expect_identical(names(flagged$sims), names(direct$sims))
  expect_identical(flagged$failures, direct$failures)
  results <- c("est_final", "post_prob_ha", "trial_success")
  expect_identical(flagged$sims[results], direct$sims[results])
})
