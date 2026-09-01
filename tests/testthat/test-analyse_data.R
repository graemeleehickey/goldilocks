test_that("analyse_data works with method = 'logrank'", {
  set.seed(6184)
  data <- data.frame(
    time = rexp(100, 0.02),
    event = sample(0:1, 100, replace = TRUE),
    treatment = rep(0:1, each = 50)
  )
  res <- analyse_data(
    data = data,
    cutpoints = NULL,
    end_of_study = 36,
    prior_surv = c(0.1, 0.1),
    N_mcmc = 10,
    single_arm = FALSE,
    method = "logrank",
    alternative = "greater",
    h0 = 0
  )
  expect_type(res, "list")
  expect_named(res, c("success", "effect"))
  expect_true(res$success >= 0 && res$success <= 1)
  expect_true(is.na(res$effect))
})

test_that("log-rank analysis rejects a nonzero h0 instead of ignoring it", {
  data <- data.frame(
    time = c(1, 2, 3, 4),
    event = c(1, 1, 1, 1),
    treatment = c(0, 0, 1, 1)
  )

  expect_error(
    analyse_data(
      data = data,
      cutpoints = NULL,
      end_of_study = 4,
      prior_surv = c(0.1, 0.1),
      N_mcmc = 10,
      single_arm = FALSE,
      method = "logrank",
      alternative = "two.sided",
      h0 = 0.1
    ),
    "'h0' must be 0 for log-rank analyses.*equal survival distributions"
  )
})

test_that("binary analyses reject missing and non-binary event indicators", {
  base_data <- data.frame(
    time = rep(36, 4),
    event = c(0, 1, 0, 1),
    treatment = c(0, 0, 1, 1)
  )
  args <- list(
    cutpoints = NULL,
    end_of_study = 36,
    prior_surv = c(0.1, 0.1),
    N_mcmc = 10,
    single_arm = FALSE,
    method = "riskdiff",
    alternative = "two.sided",
    h0 = 0
  )

  for (bad_event in list(c(0, 1, NA, 1), c(0, 1, 2, 1))) {
    data <- base_data
    data$event <- bad_event
    expect_error(
      do.call(analyse_data, c(list(data = data), args)),
      "requires binary event outcomes.*'data\\$event'.*only 0 and 1"
    )
  }
})

test_that("analyse_data works with method = 'cox'", {
  set.seed(3729)
  data <- data.frame(
    time = rexp(100, 0.02),
    event = sample(0:1, 100, replace = TRUE),
    treatment = rep(0:1, each = 50)
  )
  res <- analyse_data(
    data = data,
    cutpoints = NULL,
    end_of_study = 36,
    prior_surv = c(0.1, 0.1),
    N_mcmc = 10,
    single_arm = FALSE,
    method = "cox",
    alternative = "greater",
    h0 = 0
  )
  expect_type(res, "list")
  expect_named(res, c("success", "effect"))
  expect_true(res$success >= 0 && res$success <= 1)
  expect_type(res$effect, "double")
})

test_that("analyse_data errors clearly for non-estimable log-rank tests", {
  data <- data.frame(
    time = rep(10, 20),
    event = rep(0, 20),
    treatment = rep(0:1, each = 10)
  )

  expect_error(
    analyse_data(
      data = data,
      cutpoints = NULL,
      end_of_study = 36,
      prior_surv = c(0.1, 0.1),
      N_mcmc = 10,
      single_arm = FALSE,
      method = "logrank",
      alternative = "greater",
      h0 = 0
    ),
    "Log-rank analysis is non-estimable"
  )
})

test_that("analyse_data errors clearly for non-estimable Cox tests", {
  data <- data.frame(
    time = rep(10, 20),
    event = rep(0, 20),
    treatment = rep(0:1, each = 10)
  )

  expect_error(
    analyse_data(
      data = data,
      cutpoints = NULL,
      end_of_study = 36,
      prior_surv = c(0.1, 0.1),
      N_mcmc = 10,
      single_arm = FALSE,
      method = "cox",
      alternative = "greater",
      h0 = 0
    ),
    "Cox analysis is non-estimable"
  )
  for (engine in c("auto", "public")) {
    expect_error(
      cox_wald_test_checked(data, engine = engine),
      "Cox analysis is non-estimable"
    )
  }
})

test_that("guarded Cox engines match coxph estimates and standard errors", {
  set.seed(4291)
  data <- data.frame(
    time = c(rexp(50, rate = 0.06), rexp(50, rate = 0.04)),
    event = sample(0:1, 100, replace = TRUE),
    treatment = rep(0:1, each = 50)
  )

  automatic_fit <- cox_wald_test(data)
  public_fit <- cox_wald_test(data, engine = "public")
  survival_fit <- coxph(Surv(time, event) ~ treatment, data = data)

  expect_equal(automatic_fit$estimate, unname(survival_fit$coefficients[1]))
  expect_equal(
    automatic_fit$std_error,
    sqrt(unname(survival_fit$var[1, 1]))
  )
  expect_equal(public_fit, automatic_fit)

  compatibility <- coxph_fit_compatibility(refresh = TRUE)
  if (compatibility$compatible) {
    expect_equal(
      cox_wald_test(data, engine = "fast"),
      public_fit
    )
  }
})

test_that("cox_wald_test matches public coxph for tied event times", {
  data <- data.frame(
    time = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6),
    event = c(1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1),
    treatment = rep(0:1, 6)
  )

  fit <- cox_wald_test(data)
  reference <- survival::coxph(
    survival::Surv(time, event) ~ treatment,
    data = data,
    ties = "efron",
    singular.ok = FALSE,
    model = FALSE,
    x = FALSE,
    y = FALSE
  )

  expect_equal(fit$estimate, unname(stats::coef(reference)["treatment"]))
  expect_equal(
    fit$std_error,
    sqrt(unname(stats::vcov(reference)["treatment", "treatment"]))
  )
})

test_that("Cox compatibility resolver records and validates survival", {
  compatibility <- coxph_fit_compatibility(refresh = TRUE)

  expect_named(
    compatibility,
    c("compatible", "fitter", "version", "reason")
  )
  expect_type(compatibility$compatible, "logical")
  expect_length(compatibility$compatible, 1L)
  expect_match(compatibility$version, "^[0-9]+\\.[0-9]+")
  expect_type(compatibility$reason, "character")
  if (compatibility$compatible) {
    expect_type(compatibility$fitter, "closure")
  } else {
    expect_null(compatibility$fitter)
  }
})

test_that("Cox auto engine falls back when the fast path is incompatible", {
  set.seed(7701)
  data <- data.frame(
    time = rexp(100, rate = 0.04),
    event = rbinom(100, size = 1, prob = 0.65),
    treatment = rep(0:1, each = 50)
  )
  unavailable <- list(
    compatible = FALSE,
    fitter = NULL,
    version = "test",
    reason = "simulated incompatible signature"
  )

  expect_equal(
    cox_wald_test(data, compatibility = unavailable),
    cox_wald_test(data, engine = "public")
  )
  expect_error(
    cox_wald_test(
      data,
      engine = "fast",
      compatibility = unavailable
    ),
    "fast path is unavailable.*simulated incompatible signature"
  )
})

test_that("Cox auto engine rejects malformed fast-path results safely", {
  set.seed(7702)
  data <- data.frame(
    time = rexp(100, rate = 0.04),
    event = rbinom(100, size = 1, prob = 0.65),
    treatment = rep(0:1, each = 50)
  )
  malformed <- list(
    compatible = TRUE,
    fitter = function(...) {
      list(coefficients = 0, var = matrix(1))
    },
    version = "test",
    reason = "simulated malformed result"
  )

  expect_equal(
    cox_wald_test(data, compatibility = malformed),
    cox_wald_test(data, engine = "public")
  )
  expect_error(
    cox_wald_test(
      data,
      engine = "fast",
      compatibility = malformed
    ),
    "incompatible result structure"
  )
})

test_that("Cox analysis defines singular and convergence-failure behavior", {
  singular_data <- data.frame(
    time = seq_len(10),
    event = rep(1, 10),
    treatment = rep(0, 10)
  )
  for (engine in c("auto", "public")) {
    expect_error(
      cox_wald_test_checked(singular_data, engine = engine),
      "Cox analysis is non-estimable"
    )
  }

  separated_data <- data.frame(
    time = seq_len(20),
    event = c(rep(1, 10), rep(0, 10)),
    treatment = rep(0:1, each = 10)
  )
  for (engine in c("auto", "public")) {
    expect_error(
      cox_wald_test_checked(separated_data, engine = engine),
      "Cox analysis is non-estimable:.*reported:.*coefficient may be infinite"
    )
  }
})

test_that("analyse_data works with method = 'bayes-surv' (two-arm)", {
  set.seed(8415)
  data <- data.frame(
    time = rexp(100, 0.02),
    event = sample(0:1, 100, replace = TRUE),
    treatment = rep(0:1, each = 50)
  )
  res <- analyse_data(
    data = data,
    cutpoints = NULL,
    end_of_study = 36,
    prior_surv = c(0.1, 0.1),
    N_mcmc = 100,
    single_arm = FALSE,
    method = "bayes-surv",
    alternative = "greater",
    h0 = 0
  )
  expect_type(res, "list")
  expect_named(res, c("success", "effect"))
  expect_true(res$success >= 0 && res$success <= 1)
  expect_length(res$effect, 1)
})

test_that("Bayesian survival sufficient-statistics analysis is equivalent", {
  # These completed data sets cover both supported arm configurations and both
  # the exponential and piecewise-exponential models. Fixed seeds make this a
  # strict equivalence check of the complete posterior analysis, including the
  # event-probability conversion and posterior tail calculation.
  cases <- list(
    "two-arm, single interval" = list(
      data = data.frame(
        time = c(3, 8, 15, 4, 10, 18),
        event = c(1, 0, 1, 0, 1, 1),
        treatment = c(0, 0, 0, 1, 1, 1)
      ),
      cutpoints = NULL,
      single_arm = FALSE,
      h0 = 0
    ),
    "two-arm, piecewise" = list(
      data = data.frame(
        time = c(3, 8, 15, 4, 10, 18),
        event = c(1, 0, 1, 0, 1, 1),
        treatment = c(0, 0, 0, 1, 1, 1)
      ),
      cutpoints = c(5, 12),
      single_arm = FALSE,
      h0 = 0
    ),
    "single-arm, single interval" = list(
      data = data.frame(
        time = c(3, 8, 15, 19),
        event = c(1, 0, 1, 0),
        treatment = rep(1, 4)
      ),
      cutpoints = NULL,
      single_arm = TRUE,
      h0 = 0.5
    ),
    "single-arm, piecewise" = list(
      data = data.frame(
        time = c(3, 8, 15, 19),
        event = c(1, 0, 1, 0),
        treatment = rep(1, 4)
      ),
      cutpoints = c(5, 12),
      single_arm = TRUE,
      h0 = 0.5
    )
  )

  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    data_summ <- posterior_sufficient_stats(
      data = case$data,
      cutpoints = case$cutpoints,
      single_arm = case$single_arm
    )

    set.seed(8246)
    from_data <- analyse_data(
      data = case$data,
      cutpoints = case$cutpoints,
      end_of_study = 24,
      prior_surv = c(0.5, 0.25),
      N_mcmc = 100,
      single_arm = case$single_arm,
      method = "bayes-surv",
      alternative = "greater",
      h0 = case$h0
    )
    set.seed(8246)
    from_stats <- analyse_bayes_surv_sufficient_stats(
      data_summ = data_summ,
      cutpoints = case$cutpoints,
      end_of_study = 24,
      prior_surv = c(0.5, 0.25),
      N_mcmc = 100,
      single_arm = case$single_arm,
      alternative = "greater",
      h0 = case$h0
    )

    expect_identical(from_stats, from_data, info = case_name)
  }
})

test_that("trusted Bayesian survival analysis exactly matches checked analysis", {
  data <- data.frame(
    time = c(3, 8, 15, 4, 10, 18),
    event = c(1, 0, 1, 0, 1, 1),
    treatment = c(0, 0, 0, 1, 1, 1)
  )
  cutpoints <- c(5, 12)
  data_summ <- posterior_sufficient_stats(data, cutpoints, FALSE)
  prior <- normalize_gamma_prior(
    list(
      control = c(0.5, 0.25),
      treatment = c(1, 0.5)
    ),
    n_intervals = 3,
    single_arm = FALSE,
    name = "prior_surv"
  )
  args <- list(
    data_summ = data_summ,
    cutpoints = cutpoints,
    end_of_study = 24,
    prior_surv = prior,
    N_mcmc = 100,
    single_arm = FALSE,
    alternative = "less",
    h0 = 0
  )

  set.seed(8247)
  checked <- do.call(analyse_bayes_surv_sufficient_stats, args)
  checked_seed <- .Random.seed
  set.seed(8247)
  kernel <- do.call(analyse_bayes_surv_sufficient_stats_kernel, args)
  kernel_seed <- .Random.seed

  expect_identical(kernel, checked)
  expect_identical(
    attr(kernel, "mc_counts", exact = TRUE),
    attr(checked, "mc_counts", exact = TRUE)
  )
  expect_identical(kernel_seed, checked_seed)
})

test_that("analyse_data works with method = 'bayes-surv' and alternative = 'less'", {
  set.seed(5091)
  data <- data.frame(
    time = rexp(100, 0.02),
    event = sample(0:1, 100, replace = TRUE),
    treatment = rep(0:1, each = 50)
  )
  res_greater <- analyse_data(
    data = data,
    cutpoints = NULL,
    end_of_study = 36,
    prior_surv = c(0.1, 0.1),
    N_mcmc = 500,
    single_arm = FALSE,
    method = "bayes-surv",
    alternative = "greater",
    h0 = 0
  )
  res_less <- analyse_data(
    data = data,
    cutpoints = NULL,
    end_of_study = 36,
    prior_surv = c(0.1, 0.1),
    N_mcmc = 500,
    single_arm = FALSE,
    method = "bayes-surv",
    alternative = "less",
    h0 = 0
  )
  # P(effect > 0) + P(effect < 0) should be approximately 1
  # (exact only if no draws == 0)
  expect_equal(res_greater$success + res_less$success, 1, tolerance = 0.05)
})

test_that("analyse_data works with method = 'riskdiff'", {
  event <- c(rep(1, 10), rep(0, 40), rep(1, 20), rep(0, 30))
  data <- data.frame(
    time = rep(36, 100),
    event = event,
    treatment = rep(0:1, each = 50)
  )
  res <- analyse_data(
    data = data,
    cutpoints = NULL,
    end_of_study = 36,
    prior_surv = c(0.1, 0.1),
    N_mcmc = 10,
    single_arm = FALSE,
    method = "riskdiff",
    alternative = "greater",
    h0 = 0
  )
  expect_type(res, "list")
  expect_named(res, c("success", "effect"))
  expect_true(res$success >= 0 && res$success <= 1)
  expect_equal(res$effect, 0.2)
  expect_gt(res$success, 0.5)

  fit <- risk_difference_wald_test_checked(
    data = data,
    end_of_study = 36,
    alternative = "greater",
    h0 = 0
  )
  expect_equal(fit$estimate, 0.2)
  expect_equal(fit$variance, 0.4 * 0.6 / 50 + 0.2 * 0.8 / 50)
})

test_that("analyse_data risk-difference alternatives preserve direction", {
  data <- data.frame(
    time = rep(36, 100),
    event = c(rep(1, 10), rep(0, 40), rep(1, 20), rep(0, 30)),
    treatment = rep(0:1, each = 50)
  )
  args <- list(
    data = data,
    cutpoints = NULL,
    end_of_study = 36,
    prior_surv = c(0.1, 0.1),
    N_mcmc = 10,
    single_arm = FALSE,
    method = "riskdiff",
    h0 = 0.1
  )

  less <- do.call(analyse_data, c(args, alternative = "less"))
  greater <- do.call(analyse_data, c(args, alternative = "greater"))

  expect_equal(less$success + greater$success, 1)
  expect_lt(less$success, 0.5)
  expect_gt(greater$success, 0.5)
})

test_that("analyse_data method = 'riskdiff' rejects incomplete censored outcomes", {
  data <- data.frame(
    time = c(5, 36, 7, 36),
    event = c(0, 0, 1, 1),
    treatment = c(0, 0, 1, 1)
  )

  expect_error(
    analyse_data(
      data = data,
      cutpoints = NULL,
      end_of_study = 36,
      prior_surv = c(0.1, 0.1),
      N_mcmc = 10,
      single_arm = FALSE,
      method = "riskdiff",
      alternative = "two.sided",
      h0 = 0
    ),
    "requires all censored subjects"
  )
})

test_that("analyse_data works with method = 'bayes-bin' and Monte Carlo", {
  set.seed(4912)
  data <- data.frame(
    time = rep(36, 100),
    event = c(rep(1, 10), rep(0, 40), rep(1, 30), rep(0, 20)),
    treatment = rep(0:1, each = 50)
  )

  res <- analyse_data(
    data = data,
    cutpoints = NULL,
    end_of_study = 36,
    prior_surv = c(0.1, 0.1),
    N_mcmc = 5000,
    single_arm = FALSE,
    method = "bayes-bin",
    alternative = "greater",
    h0 = 0,
    prior_bin = c(1, 1),
    bin_method = "mc"
  )

  expect_type(res, "list")
  expect_named(res, c("success", "effect"))
  expect_true(res$success >= 0 && res$success <= 1)
  expect_length(res$effect, 1)
  expect_equal(res$effect, 20 / 52, tolerance = 0.02)
  expect_true(res$success > 0.9)
})

test_that("analyse_data method = 'bayes-bin' quadrature tails sum to 1", {
  data <- data.frame(
    time = rep(36, 40),
    event = c(rep(1, 8), rep(0, 12), rep(1, 13), rep(0, 7)),
    treatment = rep(0:1, each = 20)
  )
  args <- list(
    data = data,
    cutpoints = NULL,
    end_of_study = 36,
    prior_surv = c(0.1, 0.1),
    N_mcmc = 10,
    single_arm = FALSE,
    method = "bayes-bin",
    h0 = 0,
    prior_bin = c(1, 1),
    bin_method = "quadrature"
  )

  res_less <- do.call(analyse_data, c(args, alternative = "less"))
  res_greater <- do.call(analyse_data, c(args, alternative = "greater"))

  expect_equal(res_less$success + res_greater$success, 1, tolerance = 1e-8)
  expect_type(res_less$effect, "double")
})

test_that("analyse_data method = 'bayes-bin' engines agree on stable data", {
  data <- data.frame(
    time = rep(36, 400),
    event = c(rep(1, 80), rep(0, 120), rep(1, 120), rep(0, 80)),
    treatment = rep(0:1, each = 200)
  )
  args <- list(
    data = data,
    cutpoints = NULL,
    end_of_study = 36,
    prior_surv = c(0.1, 0.1),
    N_mcmc = 50000,
    single_arm = FALSE,
    method = "bayes-bin",
    alternative = "greater",
    h0 = 0,
    prior_bin = c(1, 1)
  )

  set.seed(3319)
  res_mc <- do.call(analyse_data, c(args, bin_method = "mc"))
  res_normal <- do.call(
    analyse_data,
    c(args, bin_method = "normal")
  )
  res_quadrature <- do.call(
    analyse_data,
    c(args, bin_method = "quadrature")
  )

  expect_equal(res_mc$success, res_quadrature$success, tolerance = 0.01)
  expect_equal(res_normal$success, res_quadrature$success, tolerance = 0.02)
})

test_that("analyse_data method = 'bayes-bin' works for single-arm data", {
  data <- data.frame(
    time = rep(36, 50),
    event = c(rep(1, 15), rep(0, 35)),
    treatment = rep(1, 50)
  )

  res <- analyse_data(
    data = data,
    cutpoints = NULL,
    end_of_study = 36,
    prior_surv = c(0.1, 0.1),
    N_mcmc = 1000,
    single_arm = TRUE,
    method = "bayes-bin",
    alternative = "less",
    h0 = 0.4,
    prior_bin = c(1, 1),
    bin_method = "quadrature"
  )

  expect_true(res$success > 0.9)
  expect_type(res$effect, "double")
})

test_that("analyse_data method = 'bayes-bin' rejects incomplete censored outcomes", {
  data <- data.frame(
    time = c(36, 5, 36, 36),
    event = c(1, 0, 1, 0),
    treatment = c(0, 0, 1, 1)
  )

  expect_error(
    analyse_data(
      data = data,
      cutpoints = NULL,
      end_of_study = 36,
      prior_surv = c(0.1, 0.1),
      N_mcmc = 1000,
      single_arm = FALSE,
      method = "bayes-bin",
      alternative = "greater",
      h0 = 0,
      prior_bin = c(1, 1),
      bin_method = "mc"
    ),
    "requires all censored subjects"
  )
})

test_that("analyse_data works with method = 'bayes-surv' (single-arm)", {
  set.seed(9473)
  data <- data.frame(
    time = rexp(50, 0.02),
    event = sample(0:1, 50, replace = TRUE),
    treatment = rep(1, 50)
  )
  res <- analyse_data(
    data = data,
    cutpoints = NULL,
    end_of_study = 36,
    prior_surv = c(0.1, 0.1),
    N_mcmc = 100,
    single_arm = TRUE,
    method = "bayes-surv",
    alternative = "greater",
    h0 = 0
  )
  expect_type(res, "list")
  expect_true(res$success >= 0 && res$success <= 1)
  expect_length(res$effect, 1)
})

test_that("analyse_data works with piecewise cutpoints", {
  set.seed(4318)
  data <- data.frame(
    time = rexp(100, 0.02),
    event = sample(0:1, 100, replace = TRUE),
    treatment = rep(0:1, each = 50)
  )
  res <- analyse_data(
    data = data,
    cutpoints = 12,
    end_of_study = 36,
    prior_surv = c(0.1, 0.1),
    N_mcmc = 100,
    single_arm = FALSE,
    method = "bayes-surv",
    alternative = "greater",
    h0 = 0
  )
  expect_type(res, "list")
  expect_true(res$success >= 0 && res$success <= 1)
})

test_that("analyse_data one-sided cox: less + greater sum to ~1", {
  set.seed(5813)
  data <- data.frame(
    time = rexp(100, 0.02),
    event = sample(0:1, 100, replace = TRUE),
    treatment = rep(0:1, each = 50)
  )
  args <- list(
    data = data,
    cutpoints = NULL,
    end_of_study = 36,
    prior_surv = c(0.1, 0.1),
    N_mcmc = 10,
    single_arm = FALSE,
    method = "cox",
    h0 = 0
  )
  res_less <- do.call(analyse_data, c(args, alternative = "less"))
  res_greater <- do.call(analyse_data, c(args, alternative = "greater"))
  res_two <- do.call(analyse_data, c(args, alternative = "two.sided"))
  # One-sided p-values should sum to 1
  # success = 1 - p, so (1 - success_less) + (1 - success_greater) = 1
  expect_equal(
    (1 - res_less$success) + (1 - res_greater$success),
    1,
    tolerance = 1e-10
  )
  # Two-sided success should be <= max of one-sided
  expect_true(res_two$success <= max(res_less$success, res_greater$success))
})

test_that("analyse_data one-sided logrank: less + greater sum to ~1", {
  set.seed(7249)
  data <- data.frame(
    time = rexp(100, 0.02),
    event = sample(0:1, 100, replace = TRUE),
    treatment = rep(0:1, each = 50)
  )
  args <- list(
    data = data,
    cutpoints = NULL,
    end_of_study = 36,
    prior_surv = c(0.1, 0.1),
    N_mcmc = 10,
    single_arm = FALSE,
    method = "logrank",
    h0 = 0
  )
  res_less <- do.call(analyse_data, c(args, alternative = "less"))
  res_greater <- do.call(analyse_data, c(args, alternative = "greater"))
  expect_equal(
    (1 - res_less$success) + (1 - res_greater$success),
    1,
    tolerance = 1e-10
  )
})

test_that("analyse_data one-sided cox detects beneficial treatment", {
  set.seed(8164)
  # Treatment has much lower hazard (beneficial)
  data <- data.frame(
    time = c(rexp(50, rate = 0.10), rexp(50, rate = 0.02)),
    event = rep(1L, 100),
    treatment = rep(0:1, each = 50)
  )
  res <- analyse_data(
    data = data,
    cutpoints = NULL,
    end_of_study = 36,
    prior_surv = c(0.1, 0.1),
    N_mcmc = 10,
    single_arm = FALSE,
    method = "cox",
    alternative = "less",
    h0 = 0
  )
  # "less" = treatment beneficial; should have high success
  expect_true(res$success > 0.95)
})

test_that("analyse_data one-sided cox uses h0 as the null log hazard ratio", {
  set.seed(123)
  data <- data.frame(
    time = c(rexp(200, rate = 0.10), rexp(200, rate = 0.12)),
    event = rep(1L, 400),
    treatment = rep(0:1, each = 200)
  )

  args <- list(
    data = data,
    cutpoints = NULL,
    end_of_study = 36,
    prior_surv = c(0.1, 0.1),
    N_mcmc = 10,
    single_arm = FALSE,
    method = "cox",
    alternative = "less"
  )
  res_superiority <- do.call(analyse_data, c(args, h0 = 0))
  res_noninferiority <- do.call(analyse_data, c(args, h0 = log(1.5)))

  expect_true(res_superiority$success < 0.05)
  expect_true(res_noninferiority$success > 0.95)
})

test_that("analyse_data validates h0 before Bayesian analysis", {
  data <- data.frame(
    time = rep(36, 4),
    event = c(0, 1, 0, 1),
    treatment = c(0, 0, 1, 1)
  )
  args <- list(
    data = data,
    cutpoints = NULL,
    end_of_study = 36,
    prior_surv = c(0.1, 0.1),
    N_mcmc = 10,
    single_arm = FALSE,
    method = "bayes-bin",
    alternative = "greater",
    prior_bin = c(1, 1),
    bin_method = "quadrature"
  )

  expect_error(
    do.call(analyse_data, c(args, h0 = NaN)),
    "single finite"
  )
  expect_error(
    do.call(analyse_data, c(args, h0 = -1.1)),
    "\\[-1, 1\\]"
  )
})
