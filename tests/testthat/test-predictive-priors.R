predictive_prior_args <- function() {
  list(
    data = data.frame(
      id = 1:20,
      treatment = 1L,
      enrollment = (0:19) / 100,
      time = 0.01,
      event = 0L,
      status = "pending"
    ),
    data_cut = 0.2,
    look = 1,
    N_total = 40,
    end_of_study = 1,
    single_arm = TRUE,
    method = "bayes-surv",
    alternative = "less",
    h0 = 0.5,
    prior_surv = c(1000, 1e6),
    prior_surv_final = c(0.1, 0.1),
    N_impute = 40,
    N_mcmc = 200,
    prob_ha = 0.975,
    seed = 9131
  )
}

test_that("both predictive sample sizes use the final Bayesian survival prior", {
  args <- predictive_prior_args()
  weak_analysis <- do.call(evaluate_interim, args)
  args$prior_surv_final <- c(1000, 1)
  strong_analysis <- do.call(evaluate_interim, args)

  # With virtually no predicted events, the weak analysis prior supports a
  # risk below 0.5. A final prior concentrated near unit event risk does not,
  # even though it must leave the low-risk predictive model unchanged.
  expect_equal(weak_analysis$probabilities$probability, c(1, 1))
  expect_equal(strong_analysis$probabilities$probability, c(0, 0))
  expect_identical(
    weak_analysis$diagnostics$posterior,
    strong_analysis$diagnostics$posterior
  )
  expect_identical(weak_analysis$decision$decision, "stop_expected_success")
  expect_identical(strong_analysis$decision$decision, "stop_futility")
})

test_that("outcome prediction still uses the interim prior", {
  args <- predictive_prior_args()
  low_risk <- do.call(evaluate_interim, args)
  args$prior_surv <- c(1000, 100)
  high_risk <- do.call(evaluate_interim, args)

  # Keep the final analysis fixed while changing the predicted event times.
  expect_equal(low_risk$probabilities$probability, c(1, 1))
  expect_equal(high_risk$probabilities$probability, c(0, 0))
  expect_identical(
    low_risk$metadata$design$prior_surv_final,
    high_risk$metadata$design$prior_surv_final
  )
})

test_that("observed interim evaluation retains and validates both Gamma priors", {
  args <- predictive_prior_args()
  args$single_arm <- FALSE
  args$data$treatment <- rep(0:1, 10)
  args$h0 <- 0
  args$cutpoints <- 0.5
  args$prior_surv <- list(control = c(1, 2), treatment = c(3, 4))
  args$prior_surv_final <- list(
    treatment = rbind(shape = c(5, 6), rate = c(10, 12)),
    control = c(2, 4)
  )
  out <- do.call(evaluate_interim, args)
  priors <- out$diagnostics$prior

  expect_identical(priors$stage, rep(c("interim", "final"), each = 4))
  expect_equal(priors$shape, c(1, 1, 3, 3, 2, 2, 5, 6))
  expect_equal(priors$rate, c(2, 2, 4, 4, 4, 4, 10, 12))
  expect_identical(out$metadata$prior_design, priors)
  expect_identical(
    out$metadata$design$prior_surv_final,
    normalize_gamma_prior(args$prior_surv_final, 2L, FALSE, "prior_surv_final")
  )

  args$prior_surv_final <- list(control = c(2, 4))
  expect_error(do.call(evaluate_interim, args), "prior_surv_final")
  args$prior_surv_final <- c(0.1, Inf)
  expect_error(do.call(evaluate_interim, args), "prior_surv_final")
  args$prior_surv_final <- matrix(1, nrow = 2, ncol = 3)
  expect_error(do.call(evaluate_interim, args), "prior_surv_final")
})

test_that("omitting the final interim-analysis prior preserves the shared prior", {
  args <- predictive_prior_args()
  args$prior_surv_final <- NULL
  implicit <- do.call(evaluate_interim, args)
  args$prior_surv_final <- args$prior_surv
  explicit <- do.call(evaluate_interim, args)

  # The recorded call differs because one supplied the default explicitly.
  implicit$metadata$call <- NULL
  explicit$metadata$call <- NULL
  expect_identical(implicit, explicit)
})

test_that("other completed-data methods ignore the final survival prior at interim", {
  args <- predictive_prior_args()
  args$single_arm <- FALSE
  args$data$treatment <- rep(0:1, 10)
  args$data$event[1:4] <- 1L
  args$data$status[1:4] <- "event"
  args$h0 <- 0
  args$prior_surv <- c(1, 10)

  for (method in c("bayes-bin", "logrank", "cox", "rmst", "riskdiff-fm")) {
    args$method <- method
    args$prior_surv_final <- c(0.1, 0.1)
    first <- do.call(evaluate_interim, args)
    args$prior_surv_final <- c(1000, 1)
    second <- do.call(evaluate_interim, args)
    expect_identical(first$probabilities, second$probabilities, info = method)
    expect_identical(first$trace, second$trace, info = method)
  }
})

test_that("trial simulations apply the final prior in their predictive decisions", {
  args <- list(
    hazard_treatment = 0,
    N_total = 40,
    interim_look = 20,
    end_of_study = 1,
    lambda = 1,
    prior_surv = c(1000, 1e6),
    method = "bayes-surv",
    alternative = "less",
    h0 = 0.5,
    N_impute = 40,
    N_mcmc = 200,
    N_trials = 2,
    seed = 9132,
    backend = "sequential",
    return_trace = TRUE
  )
  weak_analysis <- do.call(
    sim_trials,
    c(args, list(prior_surv_final = c(0.1, 0.1)))
  )
  strong_analysis <- do.call(
    sim_trials,
    c(args, list(prior_surv_final = c(1000, 1)))
  )

  expect_equal(nrow(weak_analysis$failures), 0)
  expect_equal(nrow(strong_analysis$failures), 0)
  expect_true(all(weak_analysis$traces$ppp_stop_now == 1))
  expect_true(all(weak_analysis$traces$ppp_success_at_max == 1))
  expect_true(all(strong_analysis$traces$ppp_stop_now == 0))
  expect_true(all(strong_analysis$traces$ppp_success_at_max == 0))
  expect_true(all(weak_analysis$sims$trial_success))
  expect_false(any(strong_analysis$sims$trial_success))
})
