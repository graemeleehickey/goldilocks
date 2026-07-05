# --- pwe_sim ---

test_that("pwe_sim returns correct number of rows", {
  set.seed(6274)
  out <- pwe_sim(n = 50, hazard = 0.01, cutpoints = 0)
  expect_equal(nrow(out), 50)
})

test_that("pwe_sim returns data frame with time and event columns", {
  set.seed(3891)
  out <- pwe_sim(n = 10, hazard = 0.05, cutpoints = 0, maxtime = 36)
  expect_s3_class(out, "data.frame")
  expect_named(out, c("time", "event"))
})

test_that("pwe_sim censors at maxtime", {
  set.seed(8472)
  out <- pwe_sim(n = 200, hazard = 0.001, cutpoints = 0, maxtime = 10)
  expect_true(all(out$time <= 10))
  # With low hazard, most should be censored
  expect_true(mean(out$event == 0) > 0.5)
})

test_that("pwe_sim works with piecewise hazard", {
  set.seed(5163)
  out <- pwe_sim(
    n = 100,
    hazard = c(0.005, 0.01),
    cutpoints = c(0, 12),
    maxtime = 36
  )
  expect_equal(nrow(out), 100)
  expect_true(all(out$time <= 36))
  expect_true(all(out$event %in% c(0, 1)))
})

test_that("pwe_sim without maxtime marks all as events", {
  set.seed(4029)
  out <- pwe_sim(n = 20, hazard = 0.1, cutpoints = 0)
  expect_true(all(out$event == 1))
})

test_that("pwe_sim errors on mismatched cutpoints and hazard", {
  expect_error(
    pwe_sim(n = 10, hazard = c(0.01, 0.02), cutpoints = 0),
    "length"
  )
})

test_that("pwe_sim errors when cutpoints don't start at 0", {
  expect_error(
    pwe_sim(n = 10, hazard = 0.01, cutpoints = 5),
    "First element"
  )
})

test_that("pwe_sim errors on unsorted cutpoints", {
  expect_error(
    pwe_sim(n = 10, hazard = c(0.01, 0.02, 0.03), cutpoints = c(0, 10, 5)),
    "increasing"
  )
})

test_that("pwe_sim errors on non-positive maxtime", {
  expect_error(
    pwe_sim(n = 10, hazard = 0.01, cutpoints = 0, maxtime = -1),
    "postive"
  )
})


# --- pwe_impute ---

test_that("pwe_impute returns data frame with correct dimensions", {
  set.seed(7746)
  out <- pwe_impute(
    time = c(3, 5, 8),
    hazard = c(0.01, 0.02),
    cutpoints = c(0, 12)
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
    cutpoints = c(0, 12)
  )
  expect_true(all(out$time >= obs_times))
})

test_that("pwe_impute respects maxtime censoring", {
  set.seed(9381)
  out <- pwe_impute(
    time = c(3, 5),
    hazard = c(0.002, 0.01),
    cutpoints = c(0, 12),
    maxtime = 36
  )
  expect_true(all(out$time <= 36))
})

test_that("pwe_impute without maxtime marks all as events", {
  set.seed(1847)
  out <- pwe_impute(
    time = c(1, 2, 3),
    hazard = c(0.01, 0.02),
    cutpoints = c(0, 12)
  )
  expect_true(all(out$event == 1))
})

test_that("pwe_impute errors on negative hazard", {
  expect_error(
    pwe_impute(time = 5, hazard = -0.01, cutpoints = 0),
    "less than 0"
  )
})

test_that("pwe_impute errors on negative time", {
  expect_error(
    pwe_impute(time = -1, hazard = 0.01, cutpoints = 0),
    "positive"
  )
})


# --- ppwe ---

test_that("ppwe returns probabilities in (0, 1)", {
  haz_matrix <- matrix(c(0.01, 0.02, 0.015, 0.025), ncol = 2)
  out <- ppwe(hazard = haz_matrix, end_of_study = 36, cutpoints = c(0, 12))
  expect_length(out, nrow(haz_matrix))
  expect_true(all(out > 0 & out < 1))
})

test_that("ppwe matches PWEALL across endpoint boundary cases", {
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
  cutpoints <- c(0, 6, 12)

  for (end_of_study in c(3, 6, 9, 12, 24)) {
    expected <- apply(
      haz_matrix,
      1,
      function(hazard) {
        PWEALL::pwe(
          t = end_of_study,
          rate = hazard,
          tchange = cutpoints
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
    ppwe(hazard = haz_matrix, end_of_study = 36, cutpoints = c(0, 12)),
    "do not match"
  )
})
