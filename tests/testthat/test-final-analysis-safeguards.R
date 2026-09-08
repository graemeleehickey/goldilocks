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
