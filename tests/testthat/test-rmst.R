rmst_fixture <- function() {
  data.frame(
    time = c(0, 1, 2, 2, 4, 6, 8, 8, 1, 2, 3, 4, 5, 7, 8, 8, 8),
    event = c(1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0),
    treatment = c(rep(0L, 8), rep(1L, 9))
  )
}

rmst_final_args <- function(data, imputed_final = FALSE) {
  list(
    data_in = data,
    cutpoints = NULL,
    prior_surv_final = c(0.1, 0.1),
    N_mcmc = 1,
    single_arm = FALSE,
    imputed_final = imputed_final,
    method = "rmst",
    N_impute = 5,
    alternative = "greater",
    h0 = -0.5,
    prior_bin = c(1, 1),
    bin_method = "mc",
    binary_imputation = "event-time",
    empty_interval = "prior",
    end_of_study = 8,
    rmst_tau = 6
  )
}

test_that("RMST agrees with hand calculations on completed outcomes", {
  time <- c(1, 2, 4, 5, 2, 3, 5, 5)
  event <- c(1, 1, 1, 0, 1, 1, 0, 0)
  arm <- rep(0:1, each = 4)
  fit <- rmst_estimate(time, event, arm, 5)
  expect_equal(fit$rmst, c(control = 3, treatment = 3.75))
  expect_equal(fit$arm_variance, c(control = 0.625, treatment = 0.421875))
  expect_equal(fit$estimate, 0.75)
  expect_equal(fit$variance, 1.046875)
  for (alternative in c("greater", "less", "two.sided")) {
    h0 <- -0.25
    result <- analyse_rmst(time, event, arm, 5, alternative, h0)
    p <- switch(
      alternative,
      greater = pnorm((0.75 - h0) / sqrt(1.046875), lower.tail = FALSE),
      less = pnorm((0.75 - h0) / sqrt(1.046875)),
      two.sided = 2 * pnorm(-abs((0.75 - h0) / sqrt(1.046875)))
    )
    expect_equal(result, list(success = 1 - p, effect = 0.75))
  }
})

test_that("RMST and Greenwood uncertainty agree with survRM2", {
  skip_if_not_installed("survRM2")
  d <- rmst_fixture()
  for (tau in c(3.5, 6, 8)) {
    ref <- survRM2::rmst2(d$time, d$event, d$treatment, tau = tau)
    actual <- rmst_estimate(d$time, d$event, d$treatment, tau)
    expect_equal(actual$estimate, unname(ref$unadjusted.result[1, "Est."]))
    expect_equal(
      actual$rmst,
      c(
        control = unname(ref$RMST.arm0$rmst["Est."]),
        treatment = unname(ref$RMST.arm1$rmst["Est."])
      )
    )
    expect_equal(
      actual$variance,
      unname(ref$RMST.arm0$rmst["se"]^2 + ref$RMST.arm1$rmst["se"]^2)
    )
    expect_equal(
      analyse_rmst(d$time, d$event, d$treatment, tau, "two.sided", 0)$success,
      1 - unname(ref$unadjusted.result[1, "p"])
    )
  }
})

test_that("RMST handles the fixed horizon without extrapolating positive tails", {
  d <- rmst_fixture()
  d$time[d$treatment == 0 & d$time == 8] <- 7
  expect_error(
    rmst_estimate(d$time, d$event, d$treatment, 8),
    "non-estimable.*control arm.*extrapolation"
  )
  # No remaining survival mass requires no extrapolation of a positive tail.
  d$event[d$treatment == 0] <- 1
  fit <- rmst_estimate(d$time, d$event, d$treatment, 8)
  expect_equal(unname(fit$rmst["control"]), mean(d$time[d$treatment == 0]))
  expect_true(is.finite(fit$std_error))
  expect_error(
    rmst_estimate(d$time, d$event, d$treatment, 9),
    "non-estimable.*treatment arm"
  )

  # Events at the restriction boundary have no area after them.
  d <- rmst_fixture()
  before <- rmst_estimate(d$time, d$event, d$treatment, 8)
  d$event[d$time == 8] <- 1 - d$event[d$time == 8]
  expect_equal(rmst_estimate(d$time, d$event, d$treatment, 8), before)
  d$time[d$time == 8] <- 12
  expect_equal(rmst_estimate(d$time, d$event, d$treatment, 8), before)
})

test_that("one event-free arm is allowed but zero total information is rejected", {
  time <- c(1, 2, 4, 5, 5, 5, 5, 5)
  event <- c(1, 1, 1, 0, 0, 0, 0, 0)
  arm <- rep(0:1, each = 4)
  fit <- rmst_estimate(time, event, arm, 5)
  expect_equal(unname(fit$arm_variance["treatment"]), 0)
  expect_gt(analyse_rmst(time, event, arm, 5, "greater", 0)$success, 0.5)
  expect_error(
    analyse_rmst(rep(5, 8), rep(0, 8), arm, 5, "greater", 0),
    "RMST analysis is non-estimable: total variance"
  )
  expect_error(
    analyse_rmst(rep(0, 8), rep(1, 8), arm, 5, "greater", 0),
    "RMST analysis is non-estimable: total variance"
  )
})

test_that("RMST is invariant to ordering and reverses with arm labels", {
  d <- rmst_fixture()
  a <- analyse_rmst(d$time, d$event, d$treatment, 6, "greater", -0.5)
  b <- analyse_rmst(d$time, d$event, 1 - d$treatment, 6, "less", 0.5)
  expect_equal(a$success, b$success)
  expect_equal(a$effect, -b$effect)
  d <- d[rev(seq_len(nrow(d))), ]
  expect_equal(
    analyse_rmst(d$time, d$event, d$treatment == 1, 6, "greater", -0.5),
    a
  )
})

test_that("RMST validates outcome data, horizons, margins, and configuration", {
  d <- rmst_fixture()
  expect_error(
    rmst_estimate(d$time, d$event, rep(0, nrow(d)), 6),
    "both treatment arms"
  )
  expect_error(
    rmst_estimate(d$time[-1], d$event, d$treatment, 6),
    "equally sized"
  )
  expect_error(
    rmst_estimate(replace(d$time, 1, NA), d$event, d$treatment, 6),
    "finite non-negative"
  )
  expect_error(
    rmst_estimate(d$time, replace(d$event, 1, 2), d$treatment, 6),
    "only 0 and 1"
  )
  args <- list(
    hazard_control = 0.1,
    hazard_treatment = 0.06,
    N_total = 40,
    end_of_study = 8,
    method = "rmst"
  )
  for (tau in list(0, -1, NA_real_, Inf, NULL, c(2, 4), "6")) {
    expect_error(
      do.call(survival_adapt, c(args, list(rmst_tau = tau))),
      "rmst_tau"
    )
  }
  expect_error(
    do.call(survival_adapt, c(args, list(rmst_tau = 9))),
    "must not exceed"
  )
  expect_error(
    do.call(survival_adapt, c(args, list(rmst_tau = 4, h0 = -5))),
    "h0.*rmst_tau"
  )
  expect_error(do.call(survival_adapt, c(args, list(h0 = Inf))), "h0")
  args$hazard_control <- NULL
  expect_error(do.call(survival_adapt, args), "two-armed")
  args$hazard_control <- 0.1
  args$imputed_final <- TRUE
  args$N_impute <- 1
  expect_error(do.call(survival_adapt, args), "at least two imputations")
  expect_error(
    do.call(sim_trials, c(args, list(N_trials = 1))),
    "at least two imputations"
  )
})

test_that("both completed-data routes use the same RMST horizon and margin", {
  d <- rmst_fixture()
  args <- list(
    cutpoints = c(3, 7),
    end_of_study = 8,
    prior_surv = c(0.1, 0.1),
    N_mcmc = 1,
    single_arm = FALSE,
    method = "rmst",
    alternative = "greater",
    h0 = -0.5,
    empty_interval = "prior",
    rmst_tau = 6
  )
  # The horizon is allowed to precede an analysis cutpoint.
  expected <- analyse_rmst(d$time, d$event, d$treatment, 6, "greater", -0.5)
  expect_identical(do.call(analyse_data, c(list(data = d), args)), expected)
  expect_identical(
    do.call(
      analyse_completed_survival,
      c(list(outcome = as.list(d), interval_widths = NULL), args)
    ),
    expected
  )
})

test_that("RMST final analysis retains observed censored participants", {
  d <- rmst_fixture()
  d$loss_to_fu <- d$event == 0 & d$time < 8
  d$subject_impute_success <- d$loss_to_fu
  actual <- do.call(analyse_final, rmst_final_args(d))
  expected <- analyse_rmst(d$time, d$event, d$treatment, 6, "greater", -0.5)
  expect_equal(actual, c(expected$success, expected$effect))
})

test_that("RMST final imputation pools effects and Greenwood variances", {
  d <- rmst_fixture()
  d$subject_impute_success <- d$event == 0 & d$time < 8
  args <- rmst_final_args(d, TRUE)
  set.seed(1401)
  actual <- do.call(analyse_final, args)
  set.seed(1401)
  hazards <- posterior(
    d,
    cutpoints = NULL,
    prior_surv = c(0.1, 0.1),
    N_mcmc = 5,
    single_arm = FALSE,
    empty_interval = "prior"
  )
  fits <- lapply(seq_len(5), function(j) {
    completed <- impute_data(
      d,
      hazards[j, , , drop = FALSE],
      8,
      NULL,
      "success",
      FALSE,
      "event-time"
    )
    # Independently use completed truncated times and empirical Greenwood variance.
    x <- split(pmin(completed$time, 6), completed$treatment)
    list(
      estimate = mean(x[[2]]) - mean(x[[1]]),
      variance = sum(vapply(
        x,
        function(v) sum((v - mean(v))^2) / length(v)^2,
        numeric(1)
      ))
    )
  })
  estimates <- vapply(fits, `[[`, numeric(1), "estimate")
  variances <- vapply(fits, `[[`, numeric(1), "variance")
  pooled <- pool_rubin_scalar(estimates, variances, "greater", -0.5)
  expect_equal(actual, c(pooled$success, pooled$estimate))

  d$time <- 8
  d$event <- 0
  d$subject_impute_success <- FALSE
  expect_error(
    do.call(analyse_final, rmst_final_args(d, TRUE)),
    "non-estimable: total variance"
  )
})

test_that("RMST simulations reproduce the completed-data result and saved design", {
  args <- list(
    hazard_control = 0.1,
    hazard_treatment = 0.06,
    N_total = 80,
    end_of_study = 12,
    prop_loss = 0.1
  )
  set.seed(1402)
  d <- do.call(sim_comp_data, args)
  set.seed(1402)
  actual <- do.call(
    survival_adapt,
    c(args, list(method = "rmst", rmst_tau = 8))
  )
  expected <- analyse_rmst(d$time, d$event, d$treatment, 8, "greater", 0)
  expect_equal(actual$est_final, expected$effect)
  expect_equal(actual$post_prob_ha, expected$success)
  expect_equal(attr(actual, "arguments")$rmst_tau, 8)
  expect_equal(attr(actual, "enrollment_design")$end_of_study, 12)
  set.seed(1403)
  default <- do.call(survival_adapt, c(args, list(method = "rmst")))
  set.seed(1403)
  explicit <- do.call(
    survival_adapt,
    c(args, list(method = "rmst", rmst_tau = 12))
  )
  expect_identical(default, explicit)
})

test_that("RMST interim predictions use the fixed shorter horizon in both completions", {
  data <- data.frame(
    time = c(1, 3, 8, 8, 2, 4, 0, 0),
    event = c(1, 1, 0, 0, 0, 0, 0, 0),
    treatment = rep(0:1, 4),
    subject_enrolled = c(rep(TRUE, 6), FALSE, FALSE),
    subject_impute_success = c(rep(FALSE, 4), TRUE, TRUE, FALSE, FALSE),
    subject_impute_futility = c(rep(FALSE, 6), TRUE, TRUE)
  )
  imputations <- list(
    n_draws = 1L,
    current = list(
      rows = 5:6,
      time = matrix(c(5, 8), nrow = 1),
      event = matrix(c(1L, 0L), nrow = 1)
    ),
    future = list(
      rows = 7:8,
      time = matrix(c(2, 8), nrow = 1),
      event = matrix(c(1L, 0L), nrow = 1)
    )
  )
  prepared <- prepare_predictive_outcomes(data, imputations, TRUE)
  result <- analyse_predictive_survival(
    prepared,
    imputations,
    1,
    8,
    NULL,
    NULL,
    FALSE,
    c(0.1, 0.1),
    1,
    "rmst",
    "greater",
    -0.5,
    "prior",
    TRUE,
    0.975,
    0.95,
    rmst_tau = 6
  )
  for (maximum in c(FALSE, TRUE)) {
    completed <- complete_predictive_data(
      data,
      imputations,
      1,
      include_future = maximum
    )
    if (!maximum) {
      completed <- completed[completed$subject_enrolled, ]
    }
    expected <- analyse_rmst(
      completed$time,
      completed$event,
      completed$treatment,
      6,
      "greater",
      -0.5
    )
    name <- if (maximum) "success_max" else "success_now"
    expect_identical(result[[name]], expected)
  }
})

test_that("observed RMST predictions preserve units, seeds, and event-time imputation", {
  d <- data.frame(
    id = 1:8,
    treatment = rep(0:1, 4),
    enrollment = 0:7,
    time = c(1, 2, 6, 5, 4, 3, 2, 1),
    event = c(1, 1, rep(0, 6)),
    status = c("event", "event", "censored", rep("pending", 5))
  )
  args <- list(
    data = d,
    data_cut = 8,
    look = 1,
    N_total = 12,
    end_of_study = 10,
    rmst_tau = 6,
    h0 = -0.5,
    method = "rmst",
    N_impute = 30,
    seed = 1404
  )
  set.seed(1405)
  rng <- .Random.seed
  original <- do.call(evaluate_interim, args)
  expect_identical(.Random.seed, rng)
  expect_equal(original$metadata$design$rmst_tau, 6)
  expect_identical(
    do.call(
      evaluate_interim,
      c(args, list(binary_imputation = "bernoulli"))
    )$probabilities,
    original$probabilities
  )
  scale <- 30
  args$data$time <- args$data$time * scale
  args$data$enrollment <- args$data$enrollment * scale
  args$data_cut <- args$data_cut * scale
  args$end_of_study <- args$end_of_study * scale
  args$rmst_tau <- args$rmst_tau * scale
  args$h0 <- args$h0 * scale
  args$prior_surv <- c(0.1, 0.1 * scale)
  scaled <- do.call(evaluate_interim, args)
  expect_identical(scaled$probabilities, original$probabilities)
  expect_identical(scaled$decision$decision, original$decision$decision)
  fit <- rmst_estimate(d$time, d$event, d$treatment, 1)
  rescaled <- rmst_estimate(d$time * scale, d$event, d$treatment, scale)
  expect_equal(rescaled$estimate, fit$estimate * scale)
  expect_equal(rescaled$std_error, fit$std_error * scale)
  expect_error(
    do.call(evaluate_interim, modifyList(args, list(rmst_tau = 400))),
    "must not exceed"
  )
})

test_that("RMST simulations agree across serial and PSOCK workers", {
  skip_on_cran()
  args <- list(
    hazard_control = 0.12,
    hazard_treatment = 0.07,
    N_total = 60,
    end_of_study = 12,
    rmst_tau = 8,
    method = "rmst",
    interim_look = 30,
    N_impute = 10,
    N_trials = 3,
    seed = 1406,
    imputed_final = TRUE,
    prop_loss = 0.1,
    return_trace = TRUE
  )
  serial <- do.call(sim_trials, c(args, list(backend = "sequential")))
  workers <- do.call(sim_trials, c(args, list(backend = "psock", ncores = 2)))
  expect_equal(nrow(serial$failures), 0)
  expect_identical(workers$sims, serial$sims)
  expect_identical(workers$traces, serial$traces)
  expect_identical(workers$failures, serial$failures)
  expect_equal(attr(serial, "arguments")$rmst_tau, 8)
})


test_that("completed RMST calculation agrees with Kaplan-Meier integration", {
  set.seed(1408)
  for (tau in c(1, 4, 8)) {
    for (n in c(8, 61, 200)) {
      time <- pmin(round(rexp(n, 0.2), 1), 8)
      event <- as.integer(time < 8)
      arm <- rep(0:1, length.out = n)
      auto <- rmst_estimate(time, event, arm, tau)
      km <- rmst_estimate(time, event, arm, tau, engine = "survfit")
      expect_equal(auto, km, tolerance = 1e-12)
    }
  }
  # A censoring time just below tau must use the observed-censoring route.
  d <- rmst_fixture()
  expect_equal(
    rmst_estimate(d$time, d$event, d$treatment, 8),
    rmst_estimate(d$time, d$event, d$treatment, 8, engine = "survfit")
  )
})
