# Fixed score fixtures make boundary, ordering, and counterfactual expectations
# independent of model fitting or the simulation seed.
rule_fixture <- function() {
  paths <- do.call(sim_trial_paths, path_test_args(N_trials = 1))
  for (i in 1:2) {
    paths$trials[[1]]$interims[[i]]$scores <- list(
      current = c(0.90, 0.95, 0.99, 1),
      maximum = c(0.80, 0.90, 0.95, 0.99)
    )
  }
  paths$trials[[1]]$finals$post_prob_ha <- c(0.94, 0.98, 0.99)
  paths
}

test_that("all comparisons are strict and decisions respect precedence", {
  paths <- rule_fixture()
  # At prob_ha=.95, current probability=.5 and maximum probability=.25.
  continue <- apply_stopping_rules(
    paths,
    Sn = .5,
    Qn = .5,
    Fn = .25,
    prob_ha = .95,
    return_trace = TRUE
  )
  expect_equal(continue$sims$N_enrolled, 80)
  expect_true(continue$sims$trial_success)
  expect_equal(continue$traces$ppp_stop_now, c(.5, .5))
  expect_equal(continue$traces$ppp_success_at_max, c(.25, .25))
  expect_true(all(continue$traces$decision == "continue"))
  immediate <- apply_stopping_rules(
    paths,
    Sn = .4,
    Qn = .4,
    Fn = 1,
    prob_ha = .95
  )
  expect_equal(immediate$sims$stopping_reason, "immediate_success")
  expect_true(immediate$sims$trial_success)
  expect_true(is.na(immediate$sims$post_prob_ha))
  expected <- apply_stopping_rules(
    paths,
    Sn = .4,
    Qn = .5,
    Fn = 1,
    prob_ha = .95
  )
  expect_equal(expected$sims$stopping_reason, "expected_success")
  expect_false(expected$sims$trial_success) # Current cohort final=.94, not max=.99.
  futile <- apply_stopping_rules(
    paths,
    Sn = .5,
    Qn = 1,
    Fn = .26,
    prob_ha = .95
  )
  expect_equal(futile$sims$stopping_reason, "futility")
  expect_false(futile$sims$trial_success)
  expect_true(futile$sims$futility_and_success_at_max)
  expect_equal(
    futile$sims$decision_time,
    paths$trials[[1]]$interims[[1]]$context$calendar_time
  )
})

test_that("prob_ha changes interim probabilities as well as final success", {
  paths <- rule_fixture()
  lower <- apply_stopping_rules(
    paths,
    Sn = .5,
    Fn = 0,
    prob_ha = .9,
    return_trace = TRUE
  )
  higher <- apply_stopping_rules(
    paths,
    Sn = .5,
    Fn = 0,
    prob_ha = .95,
    return_trace = TRUE
  )
  expect_equal(lower$traces$ppp_stop_now, .75)
  expect_equal(lower$sims$N_enrolled, 40)
  expect_equal(higher$sims$N_enrolled, 80)
  equal_final <- apply_stopping_rules(paths, Fn = 0, Sn = 1, prob_ha = .99)
  expect_false(equal_final$sims$trial_success)
  expect_false(equal_final$sims$success_at_max)
  at_one <- apply_stopping_rules(paths, Fn = 0, Sn = 0, prob_ha = 1)
  expect_equal(at_one$sims$N_enrolled, 80)
  expect_false(at_one$sims$trial_success)
})

test_that("look-specific thresholds stop at the first eligible look", {
  paths <- rule_fixture()
  out <- apply_stopping_rules(
    paths,
    Sn = c(1, .4),
    Fn = c(0, 0),
    prob_ha = .95,
    return_trace = TRUE
  )
  expect_equal(out$sims$N_enrolled, 60)
  expect_true(out$sims$trial_success)
  expect_equal(out$traces$decision, c("continue", "stop_expected_success"))
  expect_error(apply_stopping_rules(paths, Sn = c(.1, .2, .3)), "too many")
  expect_error(
    apply_stopping_rules(paths, Sn = .9, Qn = .8),
    "less than or equal"
  )
  expect_error(apply_stopping_rules(paths, prob_ha = NA), "prob_ha")
  expect_equal(
    apply_stopping_rules(paths, Fn = NULL)$sims,
    apply_stopping_rules(paths, Fn = 0)$sims
  )
})

test_that("replay is pure and does not copy or mutate its input bank", {
  paths <- rule_fixture()
  before <- serialize(paths, NULL)
  set.seed(515)
  state <- .Random.seed
  local_mocked_bindings(
    sim_comp_data = function(...) stop("must not simulate"),
    analyse_final = function(...) stop("must not fit"),
    evaluate_interim_evidence = function(...) stop("must not predict")
  )
  first <- apply_stopping_rules(
    paths,
    Sn = .4,
    prob_ha = .95,
    return_trace = TRUE
  )
  apply_stopping_rules(paths, Sn = 1, prob_ha = .99)
  second <- apply_stopping_rules(
    paths,
    Sn = .4,
    prob_ha = .95,
    return_trace = TRUE
  )
  expect_identical(first, second)
  expect_identical(.Random.seed, state)
  expect_identical(serialize(paths, NULL), before)
  roundtrip <- unserialize(before)
  expect_identical(
    apply_stopping_rules(roundtrip, Sn = .4, prob_ha = .95)$sims,
    first$sims
  )
  via_do_call <- do.call(apply_stopping_rules, list(paths = paths))
  expect_identical(via_do_call$call$paths, quote(paths))
})

test_that("later hypothetical failures cannot erase earlier decisions", {
  paths <- rule_fixture()
  paths$trials[[1]]$interims[2] <- list(NULL)
  paths$trials[[1]]$finals$post_prob_ha[3] <- NA_real_
  early <- apply_stopping_rules(paths, Sn = .4, Fn = 0, prob_ha = .95)
  expect_false(early$sims$trial_success)
  expect_false(early$sims$futility_and_success_at_max)
  expect_true(is.na(early$sims$success_at_max))
  expect_error(
    apply_stopping_rules(paths, Sn = 1, Fn = 0),
    "Interim calculation failed at look 2"
  )
  futile <- apply_stopping_rules(paths, Sn = 1, Fn = 1, prob_ha = .95)
  expect_false(futile$sims$trial_success)
  expect_true(is.na(futile$sims$futility_and_success_at_max))
})

test_that("candidate failures retain original trial IDs and correct denominators", {
  paths <- rule_fixture()
  second <- paths$trials[[1]]
  second$trial <- 2L
  paths$trials[[2]] <- second
  attr(paths, "arguments")$N_trials <- 2L
  attr(paths, "parallel_metadata")$tasks <- 2L
  paths$trials[[1]]$interims[1] <- list(NULL)
  expect_warning(
    out <- apply_stopping_rules(paths, Fn = 0, Sn = 1),
    "1 of 2 trials"
  )
  expect_identical(out$sims$trial, 2L)
  expect_identical(out$failures$trial, 1L)
  summary <- summarise_sims(out)
  expect_equal(summary$n_requested, 2)
  expect_equal(summary$n_failed, 1)
  expect_equal(summary$n_used, 1)
  expect_equal(summary$failure_rate, .5)
  expect_error(
    apply_stopping_rules(structure(paths, class = "list")),
    "supported"
  )
  paths$schema_version <- 99L
  expect_error(apply_stopping_rules(paths), "schema version")
})

test_that("futility measures and their denominators remain distinct", {
  paths <- rule_fixture()
  out <- apply_stopping_rules(paths, Fn = 1, Sn = 1, prob_ha = .95)
  summary <- summarise_sims(out)
  expect_equal(summary$stop_futility, 1)
  expect_equal(summary$futility_and_success_at_max, 1)
  expect_equal(summary$futility_and_success_at_max_n_used, 1)
  expect_equal(summary$futility_and_success_at_max_n_missing, 0)
  expect_lt(summary$futility_and_success_at_max_mc_lower, 1)
  paths$trials[[1]]$finals$post_prob_ha[3] <- .5
  summary <- summarise_sims(apply_stopping_rules(
    paths,
    Fn = 1,
    Sn = 1,
    prob_ha = .95
  ))
  expect_equal(summary$stop_futility, 1)
  expect_equal(summary$futility_and_success_at_max, 0)
  paths$trials[[1]]$finals$post_prob_ha[3] <- NA_real_
  summary <- summarise_sims(apply_stopping_rules(
    paths,
    Fn = 1,
    Sn = 1,
    prob_ha = .95
  ))
  expect_true(is.na(summary$futility_and_success_at_max))
  expect_equal(summary$futility_and_success_at_max_n_used, 0)
  expect_equal(summary$futility_and_success_at_max_n_missing, 1)
  expect_true(is.na(summary$futility_and_success_at_max_mc_lower))
})

test_that("posterior count diagnostics reclassify at each candidate threshold", {
  paths <- do.call(
    sim_trial_paths,
    path_test_args(method = "bayes-bin", bin_method = "mc", N_trials = 1)
  )
  evidence <- paths$trials[[1]]$interims[[1]]
  for (threshold in c(0, .5, .9, .95, 1)) {
    current <- evidence$scores$current
    counts <- attr(current, "mc_counts")
    expected <- vapply(
      seq_along(current),
      function(j) {
        analysis <- set_analysis_mc_counts(
          list(success = current[j]),
          counts$successes[j],
          counts$draws
        )
        classify_completed_analysis(analysis, threshold, .95)$uncertain
      },
      logical(1)
    )
    summary <- summarise_interim_evidence(
      evidence,
      0,
      1,
      1,
      threshold,
      .95,
      TRUE
    )
    expect_equal(summary$trace$inner_mc_uncertain_stop_now, sum(expected))
  }
})
