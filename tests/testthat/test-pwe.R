# --- pwe_sim ---

test_that("pwe_sim returns correct number of rows", {
  set.seed(6274)
  out <- pwe_sim(n = 50, hazard = 0.01, cutpoints = NULL)
  expect_equal(nrow(out), 50)
})

test_that("pwe_sim returns data frame with time and event columns", {
  set.seed(3891)
  out <- pwe_sim(n = 10, hazard = 0.05, cutpoints = NULL, maxtime = 36)
  expect_s3_class(out, "data.frame")
  expect_named(out, c("time", "event"))
})

test_that("pwe_sim censors at maxtime", {
  set.seed(8472)
  out <- pwe_sim(n = 200, hazard = 0.001, cutpoints = NULL, maxtime = 10)
  expect_true(all(out$time <= 10))
  # With low hazard, most should be censored
  expect_true(mean(out$event == 0) > 0.5)
})

test_that("pwe_sim works with piecewise hazard", {
  set.seed(5163)
  out <- pwe_sim(
    n = 100,
    hazard = c(0.005, 0.01),
    cutpoints = 12,
    maxtime = 36
  )
  expect_equal(nrow(out), 100)
  expect_true(all(out$time <= 36))
  expect_true(all(out$event %in% c(0, 1)))
})

test_that("pwe_sim without maxtime marks all as events", {
  set.seed(4029)
  out <- pwe_sim(n = 20, hazard = 0.1, cutpoints = NULL)
  expect_true(all(out$event == 1))
})

test_that("each cutpoint requires one additional hazard rate", {
  expect_no_error(pwe_sim(n = 2, hazard = 0.01))
  expect_no_error(pwe_sim(n = 2, hazard = c(0.01, 0.02), cutpoints = 5))
  expect_no_error(
    pwe_sim(n = 2, hazard = c(0.01, 0.02, 0.03), cutpoints = c(5, 10))
  )

  expect_error(
    pwe_sim(n = 10, hazard = c(0.01, 0.02), cutpoints = NULL),
    "one greater"
  )
  expect_error(
    pwe_sim(n = 10, hazard = 0.01, cutpoints = 5),
    "one greater"
  )
})

test_that("pwe_sim errors on non-positive cutpoints", {
  expect_error(
    pwe_sim(n = 10, hazard = c(0.01, 0.02), cutpoints = 0),
    "finite positive"
  )
  expect_error(
    pwe_sim(n = 10, hazard = c(0.01, 0.02), cutpoints = -1),
    "finite positive"
  )
})

test_that("pwe_sim errors on unsorted cutpoints", {
  expect_error(
    pwe_sim(n = 10, hazard = c(0.01, 0.02, 0.03), cutpoints = c(10, 5)),
    "strictly increasing"
  )
})

test_that("pwe_sim errors on non-positive maxtime", {
  expect_error(
    pwe_sim(n = 10, hazard = 0.01, cutpoints = NULL, maxtime = -1),
    "positive"
  )
})

test_that("pwe utilities require finite scalar maxtime values", {
  expect_error(
    pwe_sim(n = 10, hazard = 0.01, cutpoints = NULL, maxtime = Inf),
    "single finite positive"
  )
  expect_error(
    pwe_impute(time = 1, hazard = 0.01, cutpoints = NULL, maxtime = c(2, 3)),
    "single finite positive"
  )
  expect_error(
    pwe_sim(n = 2, hazard = c(0.01, 0.02), cutpoints = 5, maxtime = 5),
    "greater than the last cutpoint"
  )
})


# --- pwe_impute ---

test_that("pwe_impute returns data frame with correct dimensions", {
  set.seed(7746)
  out <- pwe_impute(
    time = c(3, 5, 8),
    hazard = c(0.01, 0.02),
    cutpoints = 12
  )
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 3)
  expect_named(out, c("time", "event"))
})

test_that("pwe_impute generates times after observed times", {
  set.seed(2615)
  obs_times <- c(1, 5, 10)
  out <- pwe_impute(
    time = obs_times,
    hazard = c(0.05, 0.10),
    cutpoints = 12
  )
  expect_true(all(out$time >= obs_times))
})

test_that("pwe_impute respects maxtime censoring", {
  set.seed(9381)
  out <- pwe_impute(
    time = c(3, 5),
    hazard = c(0.002, 0.01),
    cutpoints = 12,
    maxtime = 36
  )
  expect_true(all(out$time <= 36))
})

test_that("pwe_impute rejects maxtime before observed time", {
  expect_error(
    pwe_impute(time = 10, hazard = 0.1, cutpoints = NULL, maxtime = 5),
    "greater than or equal"
  )
})

test_that("pwe_impute permits maxtime equal to observed time", {
  out <- pwe_impute(time = 10, hazard = 0.1, cutpoints = NULL, maxtime = 10)

  expect_equal(out$time, 10)
  expect_equal(out$event, 0)
})

test_that("pwe_impute without maxtime marks all as events", {
  set.seed(1847)
  out <- pwe_impute(
    time = c(1, 2, 3),
    hazard = c(0.01, 0.02),
    cutpoints = 12
  )
  expect_true(all(out$event == 1))
})

test_that("pwe utilities validate finite non-negative hazard rates", {
  expect_error(
    pwe_impute(time = 5, hazard = -0.01, cutpoints = NULL),
    "finite non-negative"
  )
  expect_error(
    pwe_sim(n = 2, hazard = Inf, cutpoints = NULL, maxtime = 5),
    "finite non-negative"
  )
  expect_error(
    pwe_impute(time = 5, hazard = NaN, cutpoints = NULL),
    "finite non-negative"
  )
  expect_error(
    ppwe(matrix(-0.01, ncol = 1), end_of_study = 5, cutpoints = NULL),
    "finite non-negative"
  )
})

test_that("pwe utilities require maxtime for a zero final hazard", {
  expect_error(
    pwe_sim(n = 2, hazard = 0, cutpoints = NULL),
    "final hazard rate is zero"
  )
  expect_error(
    pwe_impute(time = 1, hazard = 0, cutpoints = NULL),
    "final hazard rate is zero"
  )

  out <- pwe_sim(n = 2, hazard = 0, cutpoints = NULL, maxtime = 5)
  expect_equal(out$time, c(5, 5))
  expect_equal(out$event, c(0, 0))

  imp <- pwe_impute(time = c(1, 2), hazard = 0, cutpoints = NULL, maxtime = 5)
  expect_equal(imp$time, c(5, 5))
  expect_equal(imp$event, c(0, 0))

  expect_no_error(
    pwe_sim(n = 2, hazard = c(0, 0.1), cutpoints = 1)
  )
  expect_no_error(
    pwe_sim(n = 2, hazard = c(0.1, 0), cutpoints = 1, maxtime = 5)
  )
})

test_that("pwe_impute errors on negative time", {
  expect_error(
    pwe_impute(time = -1, hazard = 0.01, cutpoints = NULL),
    "positive"
  )
})

test_that("conditional PWE event probabilities match the survival formula", {
  hazard <- c(0.02, 0.08, 0.01)
  cutpoints <- c(5, 12)
  time <- c(0, 1, 7, 15, 24)
  end_of_study <- 24
  interval_starts <- c(0, cutpoints)

  survival_at_time <- 1 -
    PWEALL::pwe(
      t = time,
      rate = hazard,
      tchange = interval_starts
    )$dist
  survival_at_endpoint <- 1 -
    PWEALL::pwe(
      t = end_of_study,
      rate = hazard,
      tchange = interval_starts
    )$dist
  expected <- (survival_at_time - survival_at_endpoint) / survival_at_time

  expect_equal(
    goldilocks:::pwe_conditional_event_probability(
      time = time,
      hazard = hazard,
      end_of_study = end_of_study,
      cutpoints = cutpoints
    ),
    expected,
    tolerance = 1e-12
  )
})

test_that("conditional PWE event probabilities remain stable in the tail", {
  out <- goldilocks:::pwe_conditional_event_probability(
    time = 40,
    hazard = 1,
    end_of_study = 41
  )

  expect_equal(out, 1 - exp(-1), tolerance = 1e-12)
})


# --- ppwe ---

test_that("ppwe returns probabilities in (0, 1)", {
  haz_matrix <- matrix(c(0.01, 0.02, 0.015, 0.025), ncol = 2)
  out <- ppwe(hazard = haz_matrix, end_of_study = 36, cutpoints = 12)
  expect_length(out, nrow(haz_matrix))
  expect_true(all(out > 0 & out < 1))
})

test_that("ppwe matches PWEALL after the final cutpoint", {
  haz_matrix <- matrix(
    c(
      0.01,
      0.02,
      0.03,
      0.02,
      0.03,
      0.04
    ),
    ncol = 3,
    byrow = TRUE
  )
  cutpoints <- c(6, 12)

  for (end_of_study in c(12.5, 24, 36)) {
    expected <- apply(
      haz_matrix,
      1,
      function(hazard) {
        PWEALL::pwe(
          t = end_of_study,
          rate = hazard,
          tchange = c(0, cutpoints)
        )$dist
      }
    )

    expect_equal(
      ppwe(
        hazard = haz_matrix,
        end_of_study = end_of_study,
        cutpoints = cutpoints
      ),
      expected,
      tolerance = 1e-12
    )
  }
})

test_that("ppwe errors on mismatched hazard columns and cutpoints", {
  haz_matrix <- matrix(c(0.01, 0.02), ncol = 1)
  expect_error(
    ppwe(hazard = haz_matrix, end_of_study = 36, cutpoints = 12),
    "one greater"
  )
})

test_that("piecewise helpers validate cutpoints and endpoint time", {
  expect_error(
    pwe_sim(n = 2, hazard = 0.01, cutpoints = NA_real_),
    "finite positive"
  )
  expect_error(
    pwe_sim(n = 2, hazard = c(0.01, 0.02, 0.03), cutpoints = c(6, 6)),
    "strictly increasing"
  )
  expect_error(
    ppwe(matrix(0.01, ncol = 1), end_of_study = Inf, cutpoints = NULL),
    "end_of_study"
  )
  expect_error(
    ppwe(matrix(c(0.01, 0.02), ncol = 2), end_of_study = 0, cutpoints = 12),
    "end_of_study"
  )
  expect_error(
    ppwe(matrix(c(0.01, 0.02), ncol = 2), end_of_study = 12, cutpoints = 12),
    "greater than the last cutpoint"
  )
})
