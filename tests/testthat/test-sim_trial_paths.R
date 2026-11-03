test_that("complete banks retain all looks and cohort finals", {
  paths <- do.call(sim_trial_paths, path_test_args())
  expect_s3_class(paths, "goldilocks_paths")
  expect_equal(length(paths$trials), 2L)
  expect_equal(nrow(paths$failures), 0L)
  expect_named(attr(paths, "arguments"), names(formals(sim_trial_paths)))
  for (path in paths$trials) {
    expect_equal(
      vapply(path$interims, function(x) x$context$planned_N, numeric(1)),
      c(40, 60)
    )
    expect_identical(path$finals$planned_N, c(40L, 60L, 80L))
    expect_true(all(vapply(
      path$interims,
      function(x) {
        length(x$scores$current) == 4L &&
          length(x$scores$maximum) == 4L
      },
      logical(1)
    )))
  }
  expect_output(print(paths), "2 trials; 2 interim looks")
  regenerated <- do.call(sim_trial_paths, attr(paths, "arguments"))
  expect_identical(regenerated$trials, paths$trials)
})

test_that("replay matches direct analyses with matched streams across methods", {
  methods <- c(
    "logrank",
    "cox",
    "rmst",
    "bayes-surv",
    "bayes-bin",
    "riskdiff-wald",
    "riskdiff-fm"
  )
  for (method in methods) {
    args <- path_test_args(
      method = method,
      bin_method = "mc",
      N_trials = 1,
      alternative = if (method == "rmst") "greater" else "less"
    )
    paths <- do.call(sim_trial_paths, args)
    args <- attr(paths, "arguments")
    rules <- list(
      Fn = c(0.01, 0.1),
      Sn = c(0.8, 0.5),
      Qn = c(1, 0.9),
      prob_ha = 0.9
    )
    replay <- do.call(
      apply_stopping_rules,
      c(list(paths = paths, return_trace = TRUE), rules)
    )
    reference <- direct_path_reference(args, rules)
    expect_identical(
      replay$sims$N_enrolled,
      as.integer(reference$N),
      info = method
    )
    expect_identical(
      replay$sims$trial_success,
      reference$success,
      info = method
    )
    expect_equal(
      replay$sims$post_prob_ha,
      reference$score,
      tolerance = 0,
      info = method
    )
    expect_equal(
      replay$sims$est_final,
      reference$effect,
      tolerance = 0,
      info = method
    )
    expect_equal(
      replay$traces[setdiff(names(replay$traces), "trial")],
      reference$trace,
      tolerance = 0,
      info = method
    )
    expect_equal(
      as.list(replay$sims[names(reference$calendar)]),
      reference$calendar,
      tolerance = 0,
      info = method
    )
  }
})

test_that("final imputation and nonstandard designs share the same reference", {
  cases <- list(
    list(
      method = "bayes-bin",
      bin_method = "normal",
      prop_loss = 0.1,
      imputed_final = TRUE
    ),
    list(
      method = "bayes-bin",
      bin_method = "mc",
      hazard_control = NULL,
      h0 = 0.6,
      prop_loss = 0.1,
      imputed_final = TRUE
    ),
    list(
      method = "bayes-surv",
      hazard_control = NULL,
      h0 = 0.6,
      prop_loss = 0.1,
      imputed_final = TRUE
    ),
    list(
      method = "bayes-surv",
      cutpoints = 6,
      hazard_treatment = c(0.05, 0.1),
      hazard_control = c(0.08, 0.09),
      prop_loss = c(control = 0.1, treatment = 0.2),
      imputed_final = TRUE
    ),
    list(method = "cox", prop_loss = 0.1, imputed_final = TRUE),
    list(
      method = "rmst",
      rmst_tau = 6,
      alternative = "greater",
      prop_loss = 0.1,
      imputed_final = TRUE
    ),
    list(
      method = "riskdiff-wald",
      prop_loss = 0.1,
      imputed_final = TRUE,
      rand_ratio = c(control = 1, treatment = 2),
      block = 6
    ),
    list(method = "riskdiff-fm", imputed_final = TRUE, h0 = 0.05)
  )
  for (extra in cases) {
    args <- do.call(path_test_args, c(extra, list(N_trials = 1)))
    paths <- do.call(sim_trial_paths, args)
    rules <- list(Fn = c(0, 0), Sn = c(1, 1), Qn = c(1, 1), prob_ha = 0.9)
    replay <- do.call(apply_stopping_rules, c(list(paths = paths), rules))
    reference <- direct_path_reference(attr(paths, "arguments"), rules)
    expect_equal(
      replay$sims$post_prob_ha,
      reference$score,
      tolerance = 0,
      info = extra$method
    )
    expect_equal(
      replay$sims$est_final,
      reference$effect,
      tolerance = 0,
      info = extra$method
    )
    expect_identical(
      replay$sims$trial_success,
      reference$success,
      info = extra$method
    )
  }
})

test_that("bank RNG streams are reproducible, preserved, and batch-stable", {
  set.seed(219)
  state <- .Random.seed
  args <- path_test_args()
  full <- do.call(sim_trial_paths, args)
  expect_identical(.Random.seed, state)
  first <- do.call(sim_trial_paths, path_test_args(N_trials = 1))
  second <- do.call(
    sim_trial_paths,
    path_test_args(N_trials = 1, trial_offset = 1L)
  )
  expect_identical(full$trials, c(first$trials, second$trials))
  expect_identical(.Random.seed, state)
  set.seed(219)
  sample.int(.Machine$integer.max - 1L, 1L)
  advanced <- .Random.seed
  set.seed(219)
  unseeded <- do.call(sim_trial_paths, path_test_args(seed = NULL))
  expect_identical(.Random.seed, advanced)
  regenerated <- do.call(
    sim_trial_paths,
    path_test_args(seed = attr(unseeded, "rng_metadata")$stream_seed)
  )
  expect_identical(regenerated$trials, unseeded$trials)
})

test_that("bank paths are identical across parallel backends", {
  sequential <- do.call(sim_trial_paths, path_test_args(backend = "sequential"))
  for (backend in c("psock", if (.Platform$OS.type != "windows") "fork")) {
    parallel <- do.call(
      sim_trial_paths,
      path_test_args(backend = backend, ncores = 2)
    )
    expect_identical(parallel$trials, sequential$trials, info = backend)
    expect_identical(parallel$failures, sequential$failures, info = backend)
  }
})

test_that("no-interim banks give only maximum-N final analyses", {
  paths <- do.call(sim_trial_paths, path_test_args(interim_look = NULL))
  out <- apply_stopping_rules(paths, return_trace = TRUE)
  expect_equal(out$sims$N_enrolled, rep(80, 2))
  expect_equal(nrow(out$traces), 0L)
  expect_true(all(is.na(out$sims$ppp_success)))
  expect_false(any(out$sims$futility_and_success_at_max))
})

test_that("invalid bank settings fail before simulation", {
  for (override in list(
    list(seed = -1),
    list(seed = 2^32),
    list(trial_offset = -1),
    list(N_trials = 0),
    list(N_impute = 0),
    list(mc_conf_level = 0.5),
    list(interim_look = c(40, 80)),
    list(hazard_treatment = -1),
    list(method = "logrank", imputed_final = TRUE),
    list(method = "riskdiff-fm", prop_loss = .1, imputed_final = TRUE)
  )) {
    expect_error(do.call(sim_trial_paths, do.call(path_test_args, override)))
  }
})

test_that("calculation failures do not discard later stages or unrelated scores", {
  original <- analyse_binary_counts
  local_mocked_bindings(analyse_binary_counts = function(counts, ...) {
    if (sum(counts$n_control, counts$n_treatment) == 80) {
      stop("maximum failed")
    }
    original(counts, ...)
  })
  expect_warning(
    paths <- do.call(sim_trial_paths, path_test_args()),
    "path calculations failed"
  )
  expect_true(all(vapply(
    paths$trials,
    function(path) {
      length(path$interims) == 2 &&
        all(vapply(
          path$interims,
          function(x) length(x$scores$current) == 4L,
          logical(1)
        ))
    },
    logical(1)
  )))
  expect_true(all(paths$failures$stage == "interim_maximum"))
  # Maximum final analyses use the general final-analysis path and still work.
  expect_silent(out <- apply_stopping_rules(paths, Fn = 0, Sn = 1, Qn = 1))
  expect_equal(nrow(out$sims), 2)
  expect_error(
    apply_stopping_rules(paths, Fn = 1, Sn = 1, Qn = 1),
    "All 2 trials failed"
  )
})

test_that("failed final analyses and generation are recorded without silent omission", {
  local_mocked_bindings(analyse_final = function(...) stop("final unavailable"))
  expect_warning(
    paths <- do.call(sim_trial_paths, path_test_args()),
    "path calculations failed"
  )
  expect_equal(sum(paths$failures$stage == "final"), 6L)
  expect_error(
    apply_stopping_rules(paths, Fn = 0, Sn = 1),
    "Required final analysis failed"
  )
  expect_silent(out <- apply_stopping_rules(paths, Fn = 1, Sn = 1, prob_ha = 1))
  expect_true(all(out$sims$stop_futility == 1))
  expect_true(all(is.na(out$sims$futility_and_success_at_max)))
  expect_equal(summarise_sims(out)$futility_and_success_at_max_n_missing, 2)
  local_mocked_bindings(sim_comp_data = function(...) {
    stop("generation unavailable")
  })
  expect_warning(
    paths <- do.call(sim_trial_paths, path_test_args()),
    "path calculations failed"
  )
  expect_equal(nrow(paths$failures), 2L)
  expect_true(all(paths$failures$stage == "generation"))
})
