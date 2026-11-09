horizon_test_args <- function(...) {
  args <- list(
    hazard_treatment = 0.05,
    hazard_control = 0.08,
    N_total = 80,
    interim_look = c(40, 60),
    lambda = 3,
    final_horizon = 18,
    endpoint_horizons = c(12, 18),
    alternative = "less",
    method = "bayes-bin",
    bin_method = "quadrature",
    N_impute = 4,
    N_mcmc = 30,
    N_trials = 2,
    seed = 914
  )
  extra <- list(...)
  args[names(extra)] <- extra
  args
}

horizon_test_schedule <- function() {
  data.frame(N = c(40, 60), horizon = c(12, 18))
}

test_that("endpoint recensoring preserves event and dropout ordering", {
  data <- data.frame(
    time = c(15, 15, 8, 12, 18),
    event = c(1, 0, 0, 1, 0),
    loss_to_fu = c(FALSE, TRUE, TRUE, FALSE, FALSE)
  )
  cut <- censor_horizon_data(data, 12)
  expect_equal(cut$time, c(12, 12, 8, 12, 12))
  expect_equal(cut$event, c(0, 0, 0, 1, 0))
  expect_identical(cut$loss_to_fu, c(FALSE, FALSE, TRUE, FALSE, FALSE))
  expect_identical(censor_horizon_data(data, 18), data)
})

test_that("late observed data inform futility but cannot inform earlier success", {
  bank <- do.call(sim_horizon_paths, horizon_test_args(N_trials = 1))
  args <- attr(bank, "arguments")
  config12 <- horizon_analysis_config(args, 12)
  config18 <- horizon_analysis_config(args, 18)
  # A single shared hazard interval crosses month 12: mere re-analysis of an
  # 18-month posterior would fail this test even without an extra cutpoint.
  data <- data.frame(
    id = 1:80,
    treatment = rep(c(0, 1), 40),
    time = rep(18, 80),
    event = rep(0, 80),
    enrollment = c(rep(0, 39), 20:60),
    loss_to_fu = FALSE
  )
  late <- data
  late$time[1:20] <- 15
  late$event[1:20] <- 1
  evaluate <- function(data, config) {
    set.seed(45)
    horizon_interim_evidence(data, config, 1, 40)
  }
  expect_identical(evaluate(data, config12), evaluate(late, config12))
  original18 <- evaluate(data, config18)
  changed18 <- evaluate(late, config18)
  expect_false(identical(original18$scores$maximum, changed18$scores$maximum))
  expect_equal(sum(changed18$diagnostics$posterior$observed_events), 20)
  expect_equal(sum(original18$diagnostics$posterior$observed_events), 0)
  # Outcomes after the calendar cut, including outcomes of future recruits,
  # cannot affect either predictive calculation.
  future <- late
  future$time[40:80] <- 2
  future$event[40:80] <- 1
  expect_identical(evaluate(late, config12), evaluate(future, config12))
  expect_identical(evaluate(late, config18), evaluate(future, config18))
})

test_that("each horizon uses the leading master-prior intervals", {
  args <- horizon_test_args(
    cutpoints = c(6, 12, 15),
    generation_cutpoints = NULL,
    prior_surv = matrix(c(1, 10, 2, 20, 3, 30, 4, 40), nrow = 2),
    N_trials = 1
  )
  bank <- do.call(sim_horizon_paths, args)
  config <- horizon_analysis_config(attr(bank, "arguments"), 12)
  expect_equal(config$cutpoints, 6)
  expect_equal(
    unname(config$prior_surv[,, "control"]),
    matrix(c(1, 10, 2, 20), 2)
  )
  posterior <- bank$trials[[1]]$interims[[1]][[1]]$diagnostics$posterior
  expect_setequal(posterior$interval_end, c(6, 12))
  expect_setequal(posterior$prior_shape, c(1, 2))
  args$prior_surv[, 3:4] <- 100
  other <- do.call(sim_horizon_paths, args)
  expect_identical(
    bank$trials[[1]]$interims[[1]][[1]],
    other$trials[[1]]$interims[[1]][[1]]
  )
})

test_that("direct horizon trials agree with full-bank replay for every method", {
  for (method in c(
    "logrank",
    "cox",
    "rmst",
    "bayes-surv",
    "bayes-bin",
    "riskdiff-wald",
    "riskdiff-fm"
  )) {
    args <- horizon_test_args(
      method = method,
      N_trials = 1,
      alternative = if (method == "rmst") "greater" else "less"
    )
    bank <- do.call(sim_horizon_paths, args)
    schedule <- horizon_test_schedule()
    for (rules in list(
      list(Fn = 0, Sn = 1),
      list(Fn = 1, Sn = 1, prob_ha = 1),
      list(Fn = 0, Sn = 0, Qn = 1, prob_ha = 0),
      list(Fn = 0, Sn = 0, Qn = 0, prob_ha = 0)
    )) {
      replay <- do.call(
        apply_horizon_rules,
        c(
          list(paths = bank, endpoint_schedule = schedule, return_trace = TRUE),
          rules
        )
      )
      direct_args <- args
      direct_args$interim_look <- NULL
      direct <- do.call(
        sim_horizon_trials,
        c(
          direct_args,
          list(endpoint_schedule = schedule, return_trace = TRUE),
          rules
        )
      )
      expect_identical(direct$sims, replay$sims, info = method)
      expect_identical(direct$traces, replay$traces, info = method)
    }
  }
})

test_that("constant horizons reproduce fixed-horizon paths with matched streams", {
  for (method in c("logrank", "bayes-surv", "bayes-bin", "rmst")) {
    args <- horizon_test_args(
      endpoint_horizons = 18,
      method = method,
      N_trials = 1,
      prop_loss = 0.1,
      imputed_final = method != "logrank",
      alternative = if (method == "rmst") "greater" else "less"
    )
    bank <- do.call(sim_horizon_paths, args)
    old_args <- args[setdiff(
      names(args),
      c("endpoint_horizons", "final_horizon")
    )]
    old_args$end_of_study <- 18
    fixed <- do.call(sim_trial_paths, old_args)
    for (i in 1:2) {
      expect_identical(
        bank$trials[[1]]$interims[[i]][[1]],
        fixed$trials[[1]]$interims[[i]],
        info = method
      )
    }
    expect_identical(
      bank$trials[[1]]$finals[names(fixed$trials[[1]]$finals)],
      fixed$trials[[1]]$finals,
      info = method
    )
  }
})

test_that("single-arm horizon-specific hypotheses and imputed finals are retained", {
  for (method in c("bayes-bin", "bayes-surv")) {
    settings <- data.frame(horizon = c(12, 18), h0 = c(0.4, 0.6))
    args <- horizon_test_args(
      hazard_control = NULL,
      method = method,
      horizon_settings = settings,
      prop_loss = 0.2,
      dropout_horizon = 12,
      imputed_final = TRUE,
      N_trials = 1
    )
    bank <- do.call(sim_horizon_paths, args)
    replay <- apply_horizon_rules(bank, horizon_test_schedule(), Fn = 0, Sn = 1)
    expect_equal(replay$sims$margin, 0.6)
    expect_equal(replay$sims$N_control, 0)
    expect_true(is.finite(replay$sims$post_prob_ha))
    args$interim_look <- NULL
    direct <- do.call(
      sim_horizon_trials,
      c(args, list(endpoint_schedule = horizon_test_schedule(), Fn = 0, Sn = 1))
    )
    expect_identical(direct$sims, replay$sims)
  }
})

test_that("a stop locks the endpoint and keeps counterfactual futility separate", {
  bank <- do.call(sim_horizon_paths, horizon_test_args(N_trials = 1))
  bank$trials[[1]]$interims[[1]][[1]]$scores$current[] <- 0.96
  bank$trials[[1]]$interims[[1]][[2]]$scores$maximum[] <- 0.1
  bank$trials[[1]]$finals$post_prob_ha <- c(0.1, 0.99, 0.1, 0.99, 0.1, 0.99)
  schedule <- horizon_test_schedule()
  out <- apply_horizon_rules(
    bank,
    schedule,
    Sn = 0.9,
    Fn = 1,
    return_trace = TRUE
  )
  expect_equal(out$sims$N_enrolled, 40)
  expect_equal(out$sims$endpoint_horizon, 12)
  expect_false(out$sims$trial_success)
  expect_true(out$sims$success_at_max)
  expect_equal(out$sims$stop_expected_success, 1)
  expect_equal(nrow(out$traces), 1)
  expect_equal(out$traces$success_data_horizon, 12)
  expect_equal(out$traces$futility_data_horizon, 18)
  expect_equal(
    out$sims$planned_collection_time - out$sims$planned_completion_time,
    6
  )
  expect_false(out$sims$futility_and_success_at_max)
  # Equality does not cross the success boundary; futility then applies.
  out <- apply_horizon_rules(bank, schedule, Sn = 1, Fn = 0.1)
  expect_equal(out$sims$stop_futility, 1)
  expect_true(out$sims$futility_and_success_at_max)
  expect_equal(summarise_sims(out)$futility_and_success_at_max, 1)
})

test_that("a continued trial changes endpoint before evaluating the scheduled look", {
  bank <- do.call(sim_horizon_paths, horizon_test_args(N_trials = 1))
  bank$trials[[1]]$interims[[1]][[1]]$scores$current[] <- 0.1
  bank$trials[[1]]$interims[[2]][[1]]$scores$current[] <- 0.1
  bank$trials[[1]]$interims[[2]][[2]]$scores$current[] <- 0.99
  out <- apply_horizon_rules(
    bank,
    horizon_test_schedule(),
    Fn = 0,
    return_trace = TRUE
  )
  expect_equal(out$sims$N_enrolled, 60)
  expect_equal(out$sims$endpoint_horizon, 18)
  expect_equal(out$traces$success_horizon, c(12, 18))
  never <- apply_horizon_rules(bank, horizon_test_schedule(), Fn = 0, Sn = 1)
  expect_equal(never$sims$N_enrolled, 80)
  expect_equal(
    never$sims$planned_completion_time,
    never$sims$accrual_stop_time + 18
  )
})

test_that("unused failures do not invalidate known horizon outcomes", {
  bank <- do.call(sim_horizon_paths, horizon_test_args(N_trials = 1))
  bank$trials[[1]]$interims[[1]][[1]]$scores$current[] <- 0.99
  bank$trials[[1]]$interims[[1]][[2]] <- NULL
  # Preserve the horizon slot when recording a failed model.
  bank$trials[[1]]$interims[[1]][2] <- list(NULL)
  bank$trials[[1]]$finals$post_prob_ha[6] <- NA_real_
  out <- apply_horizon_rules(bank, horizon_test_schedule(), Qn = 0.9)
  expect_true(out$sims$trial_success)
  expect_true(is.na(out$sims$success_at_max))
  expect_false(out$sims$futility_and_success_at_max)
  expect_error(
    apply_horizon_rules(bank, horizon_test_schedule(), Sn = 1, Fn = 1),
    "F.*prediction|Maximum-cohort"
  )
  bank$trials[[1]]$interims[[1]][[1]]$scores$current[] <- 0.1
  bank$trials[[1]]$interims[[1]][2] <- list(list(
    scores = list(maximum = rep(0.1, 4)),
    diagnostics = list()
  ))
  out <- apply_horizon_rules(bank, horizon_test_schedule(), Sn = 1, Fn = 1)
  expect_true(is.na(out$sims$futility_and_success_at_max))
  expect_equal(summarise_sims(out)$futility_and_success_at_max_n_missing, 1)
})

test_that("direct simulation skips unreached interim fits", {
  original <- horizon_interim_evidence
  calls <- integer()
  local_mocked_bindings(horizon_interim_evidence = function(
    data,
    config,
    look,
    n
  ) {
    calls <<- c(calls, look)
    original(data, config, look, n)
  })
  args <- horizon_test_args(N_trials = 1)
  args$interim_look <- NULL
  out <- do.call(
    sim_horizon_trials,
    c(
      args,
      list(
        endpoint_schedule = horizon_test_schedule(),
        Sn = 1,
        Fn = 1,
        prob_ha = 1
      )
    )
  )
  expect_equal(out$sims$stop_futility, 1)
  expect_identical(calls, c(1L, 1L))
})

test_that("horizon banks are reproducible, serializable, and batch-stable", {
  set.seed(23)
  before <- .Random.seed
  bank <- do.call(sim_horizon_paths, horizon_test_args())
  expect_identical(.Random.seed, before)
  expect_named(attr(bank, "arguments"), names(formals(sim_horizon_paths)))
  regenerated <- do.call(sim_horizon_paths, attr(bank, "arguments"))
  expect_identical(regenerated$trials, bank$trials)
  first <- do.call(sim_horizon_paths, horizon_test_args(N_trials = 1))
  second <- do.call(
    sim_horizon_paths,
    horizon_test_args(N_trials = 1, trial_offset = 1)
  )
  expect_identical(bank$trials, c(first$trials, second$trials))
  file <- tempfile(fileext = ".rds")
  on.exit(unlink(file))
  saveRDS(bank, file)
  restored <- readRDS(file)
  original <- apply_horizon_rules(
    bank,
    horizon_test_schedule(),
    return_trace = TRUE
  )
  replay <- apply_horizon_rules(
    restored,
    horizon_test_schedule(),
    return_trace = TRUE
  )
  expect_identical(replay$sims, original$sims)
  expect_identical(replay$traces, original$traces)
  expect_identical(.Random.seed, before)
  expect_output(print(bank), "2 trials; 2 interim looks; horizons: 12, 18")
})

test_that("horizon banks agree across computation methods", {
  sequential <- do.call(
    sim_horizon_paths,
    horizon_test_args(backend = "sequential")
  )
  for (backend in c("psock", if (.Platform$OS.type != "windows") "fork")) {
    parallel <- do.call(
      sim_horizon_paths,
      horizon_test_args(backend = backend, ncores = 2)
    )
    expect_identical(parallel$trials, sequential$trials, info = backend)
    expect_identical(parallel$failures, sequential$failures, info = backend)
  }
})

test_that("dropout probability retains its fixed reference horizon", {
  a <- horizon_test_args(prop_loss = 0.2, dropout_horizon = 12, N_trials = 1)
  b <- horizon_test_args(
    prop_loss = -expm1(log1p(-0.2) * 18 / 12),
    dropout_horizon = 18,
    N_trials = 1
  )
  expect_identical(
    do.call(sim_horizon_paths, a)$trials,
    do.call(sim_horizon_paths, b)$trials
  )
})

test_that("invalid horizon settings fail before simulation", {
  for (override in list(
    list(endpoint_horizons = c(18, 12)),
    list(endpoint_horizons = c(12, 12, 18)),
    list(final_horizon = 12),
    list(followup_horizon = 12),
    list(dropout_horizon = 0),
    list(interim_look = c(40, 80)),
    list(cutpoints = 18),
    list(horizon_settings = data.frame(horizon = 12)),
    list(horizon_settings = data.frame(horizon = c(12, 18), h0 = c(0, NA))),
    list(
      horizon_settings = data.frame(horizon = c(12, 18), rmst_tau = c(18, 18))
    ),
    list(method = "riskdiff-fm", imputed_final = TRUE, prop_loss = 0.1)
  )) {
    expect_error(do.call(
      sim_horizon_paths,
      do.call(horizon_test_args, override)
    ))
  }
  bank <- do.call(sim_horizon_paths, horizon_test_args())
  for (schedule in list(
    data.frame(N = c(40, 60), horizon = c(18, 12)),
    data.frame(N = c(40, 60), horizon = c(6, 18)),
    data.frame(N = c(40, 80), horizon = c(12, 18))
  )) {
    expect_error(apply_horizon_rules(bank, schedule), "endpoint_schedule")
  }
  malformed <- bank
  malformed$trials[[1]]$finals$horizon[1] <- 18
  expect_error(
    apply_horizon_rules(malformed, horizon_test_schedule()),
    "Inconsistent"
  )
})

test_that("designs without interims use the final endpoint only", {
  bank <- do.call(sim_horizon_paths, horizon_test_args(interim_look = NULL))
  schedule <- data.frame(N = numeric(), horizon = numeric())
  out <- apply_horizon_rules(bank, schedule, return_trace = TRUE)
  expect_equal(out$sims$endpoint_horizon, c(18, 18))
  expect_equal(nrow(out$traces), 0)
  args <- horizon_test_args(interim_look = NULL)
  args$interim_look <- NULL
  direct <- do.call(
    sim_horizon_trials,
    c(args, list(endpoint_schedule = schedule))
  )
  expect_identical(out$sims, direct$sims)
})
