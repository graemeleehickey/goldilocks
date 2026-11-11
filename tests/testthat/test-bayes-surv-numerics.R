# Fixed posterior draws isolate numerical accuracy from posterior sampling.
bayes_surv_numerical_fixture <- function(treatment, control = NULL) {
  array(
    c(treatment, if (is.null(control)) rep(NA_real_, length(treatment)) else control),
    dim = c(length(treatment), 1L, 2L)
  )
}

test_that("Bayesian survival comparisons retain ordering when risks round to one", {
  post <- bayes_surv_numerical_fixture(
    c(41, 40, 1001, 1000, 40),
    c(40, 41, 1000, 1001, 40)
  )
  for (alternative in c("less", "greater")) {
    out <- analyse_bayes_surv_posterior_draws(
      post, NULL, 1, FALSE, alternative, 0
    )
    expect_equal(out$success, 2 / 5)
    expect_identical(attr(out, "mc_counts"), list(successes = 2L, draws = 5L))
  }

  # The last two unequal risks differ by less than a representable double.
  # Their reported difference can be zero, but their ordering remains known.
  risks <- haz_to_prop(post, NULL, 1, FALSE)
  reference <- exp(-40) * (1 - exp(-1))
  expect_equal(risks$effect[1] / reference, 1, tolerance = 1e-12)
  expect_equal(risks$effect[2] / reference, -1, tolerance = 1e-12)
  expect_equal(risks$effect[3:5], c(0, 0, 0))
})

test_that("nonzero Bayesian survival margins retain their risk-difference scale", {
  post <- bayes_surv_numerical_fixture(
    c(41, 40, 1001, 1000, 40),
    c(40, 41, 1000, 1001, 40)
  )
  # The first two differences are approximately +/-2.68547e-18.
  margins <- c(-3e-18, -2e-18, 2e-18, 3e-18)
  expected_less <- c(0, 1, 4, 5) / 5
  for (i in seq_along(margins)) {
    for (alternative in c("less", "greater")) {
      out <- analyse_bayes_surv_posterior_draws(
        post, NULL, 1, FALSE, alternative, margins[i]
      )
      expected <- if (alternative == "less") expected_less[i] else 1 - expected_less[i]
      expect_equal(out$success, expected)
    }
  }
  for (margin in c(-1, 1) * .Machine$double.xmin * .Machine$double.eps) {
    out <- analyse_bayes_surv_posterior_draws(post, NULL, 1, FALSE, "less", margin)
    expect_equal(out$success, if (margin < 0) 1 / 5 else 4 / 5)
  }
})

test_that("Bayesian survival comparisons are strict at nonzero and boundary margins", {
  for (direction in c(-1, 1)) {
    post <- if (direction == 1) {
      bayes_surv_numerical_fixture(log(2), 0)
    } else {
      bayes_surv_numerical_fixture(0, log(2))
    }
    for (alternative in c("less", "greater")) {
      out <- analyse_bayes_surv_posterior_draws(
        post, NULL, 1, FALSE, alternative, direction / 2
      )
      expect_equal(out$success, 0)
    }
  }

  post <- bayes_surv_numerical_fixture(c(1000, 0, 0), c(0, 1000, 0))
  for (margin in c(-1, 1)) {
    for (alternative in c("less", "greater")) {
      out <- analyse_bayes_surv_posterior_draws(
        post, NULL, 1, FALSE, alternative, margin
      )
      expect_equal(out$success, as.numeric((margin == 1) == (alternative == "less")))
    }
  }
})

test_that("single-arm Bayesian survival uses strict probability thresholds", {
  post <- bayes_surv_numerical_fixture(c(0, log(2), 40, 1000))
  margins <- c(0, .5, 1 - .Machine$double.eps, 1)
  expected_less <- c(0, 1, 2, 4) / 4
  expected_greater <- c(3, 2, 2, 0) / 4
  for (i in seq_along(margins)) {
    for (alternative in c("less", "greater")) {
      out <- analyse_bayes_surv_posterior_draws(
        post, NULL, 1, TRUE, alternative, margins[i]
      )
      expected <- if (alternative == "less") expected_less[i] else expected_greater[i]
      expect_equal(out$success, expected)
    }
  }
})

test_that("prepared and general survival analyses use the same endpoint comparison", {
  # Widths (1, 2) give cumulative hazards (42, 62) and (81, 41).
  # Crossing interval hazards must be compared after integration to the endpoint.
  post <- array(c(40, 60, 1, 1, 1, 1, 40, 20), dim = c(2, 2, 2))
  local_mocked_bindings(posterior_from_prepared_stats = function(...) post)
  for (alternative in c("less", "greater")) {
    for (margin in c(-1e-19, 0, 1e-19)) {
      general <- analyse_bayes_surv_posterior_draws(
        post, 1, 3, FALSE, alternative, margin
      )
      prepared <- analyse_prepared_bayes_surv(
        data_summ = NULL, cutpoints = 1, end_of_study = 3,
        prior_surv = NULL, N_mcmc = 2, single_arm = FALSE,
        alternative = alternative, h0 = margin, interval_widths = c(1, 2)
      )
      expect_identical(prepared, general)
      if (margin == 0) expect_equal(prepared$success, .5)
    }
  }
})

test_that("identical Gamma posteriors retain symmetry at high endpoint risks", {
  set.seed(819)
  post <- array(rgamma(200000, shape = 21, rate = 21), c(100000, 1, 2))
  out <- analyse_bayes_surv_posterior_draws(post, NULL, 50, FALSE, "less", 0)
  expect_equal(out$success, mean(post[, 1, 1] < post[, 1, 2]))
  expect_equal(out$success, .5, tolerance = .01)
})

test_that("stable survival comparisons agree with ordinary risk calculations", {
  set.seed(820)
  post <- array(rgamma(4000, shape = 4, rate = 100), c(1000, 2, 2))
  treatment <- 1 - exp(-(5 * post[, 1, 1] + 7 * post[, 2, 1]))
  control <- 1 - exp(-(5 * post[, 1, 2] + 7 * post[, 2, 2]))
  for (single_arm in c(FALSE, TRUE)) {
    effect <- if (single_arm) treatment else treatment - control
    margins <- if (single_arm) c(.2, .5, .8) else c(-.2, 0, .2)
    for (margin in margins) {
      for (alternative in c("less", "greater")) {
        before <- .Random.seed
        out <- analyse_bayes_surv_posterior_draws(
          post, 5, 12, single_arm, alternative, margin
        )
        expected <- if (alternative == "less") mean(effect < margin) else mean(effect > margin)
        expect_equal(out$success, expected)
        expect_equal(out$effect, mean(effect), tolerance = 1e-12)
        expect_identical(.Random.seed, before)
      }
    }
  }
})

test_that("non-finite cumulative hazards cannot silently determine success", {
  post <- bayes_surv_numerical_fixture(1e308, 1e308)
  expect_error(
    analyse_bayes_surv_posterior_draws(post, NULL, 2, FALSE, "less", 0),
    "finite, non-negative cumulative hazards"
  )
})
