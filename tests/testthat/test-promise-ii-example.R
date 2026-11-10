source(system.file("examples", "promise-ii.R", package = "goldilocks"), local = TRUE)

test_that("PROMISE II hazards and reference tables are internally consistent", {
  design <- promise_ii_design()
  expect_equal(exp(-cumsum(design$base_hazard * c(30, 30, 30, 90))),
               c(.85, .70, .62, .54), tolerance = 2e-7)
  for (afs in promise_ii_reference()$afs) {
    expect_equal(1 - ppwe(matrix(promise_ii_hazard(afs), nrow = 1),
                         180, c(30, 60, 90)), afs)
  }
  ref <- promise_ii_reference()
  probabilities <- as.matrix(ref[, 4:12])
  expect_equal(rowSums(probabilities), rep(1, 8), tolerance = 2e-4)
  implied_n <- drop(probabilities %*% c(60, 60, 75, 75, 90, 90, 105, 105, 120))
  expect_lt(max(abs(implied_n - ref$mean_n)), .06)
})

test_that("PROMISE II data cuts retain events, dropout, and endpoint boundaries", {
  data <- data.frame(id = 1:6, treatment = 1L, enrollment = c(0, 2, 4, 6, 8, 100),
                     event_time = c(30, 60, 90, 180, 185, 15),
                     dropout_time = c(Inf, Inf, Inf, Inf, Inf, 10))
  interim <- promise_ii_cut(data, 6)
  expect_equal(interim$event, c(1, 1, 1, 0, 0, 0))
  expect_equal(interim$time, c(30, 60, 90, 94, 92, 0))
  final <- promise_ii_cut(data, 6, final = TRUE)
  expect_equal(final$event, c(1, 1, 1, 1, 0, 0))
  expect_equal(final$time, c(30, 60, 90, 180, 180, 10))
  expect_equal(final$status, c(rep("event", 4), "complete", "censored"))
  expect_equal(attr(final, "data_cut"), 280)
})

test_that("PROMISE II posterior updates use right-closed event intervals", {
  data <- data.frame(time = c(0, 30, 60, 90, 180, 10),
                     event = c(0, 1, 1, 1, 1, 0))
  sufficient <- promise_ii_sufficient(data)
  expect_equal(unname(sufficient$shape), rep(1.001, 4))
  expect_equal(unname(sufficient$rate), c(130, 90, 60, 90) + .001)
})

test_that("PROMISE II final imputation agrees with an analytical one-missing result", {
  data <- data.frame(time = c(rep(180, 38), rep(20, 20), 75),
                     event = c(rep(0, 38), rep(1, 20), 0))
  post <- promise_ii_sufficient(data)
  remaining <- c(0, 0, 15, 90)
  # Gamma Laplace transform integrates the missing subject's event probability.
  p <- 1 - prod((post$rate / (post$rate + remaining))^post$shape)
  scores <- promise_ii_score(c(20, 21), nrow(data))
  exact <- (1 - p) * scores[1] + p * scores[2]
  set.seed(4271)
  estimate <- promise_ii_final(data, 50000)
  se <- abs(diff(scores)) * sqrt(p * (1 - p) / 50000)
  expect_lt(abs(estimate - exact), 5 * se)
  # Complete data use the exact beta tail; final success is not mean(score > .977).
  complete <- data[-nrow(data), ]
  expect_equal(promise_ii_final(complete, 10), promise_ii_score(20, 58))
})

test_that("PROMISE II count predictions agree with the package", {
  set.seed(8)
  data <- promise_ii_generate(.70)
  for (look in c(1L, 4L)) {
    cut <- promise_ii_cut(data, promise_ii_design()$interim_look[look])
    set.seed(331 + look)
    actual <- promise_ii_predict(cut, 30000)
    reference <- promise_ii_package_interim(cut, look, 30000, seed = 121 + look)
    expected <- reference$probabilities$probability
    se <- sqrt((actual * (1 - actual) + expected * (1 - expected)) / 30000)
    expect_true(all(abs(actual - expected) < 5 * se + 1 / 30000))
  }
})

test_that("PROMISE II uses strict boundaries and keeps final inference after futility", {
  expect_equal(promise_ii_decision(c(.98, .05), 1), "continue")
  expect_equal(promise_ii_decision(c(.99, .01), 1), "promising")
  expect_equal(promise_ii_decision(c(.98, .04), 1), "futile")
  design <- promise_ii_design()
  design$Sn[] <- 1
  design$Fn[] <- 1
  design$prob_ha <- 0
  # Force PPmax < 1 without changing the final criterion, isolating the stop rule.
  controlled <- new.env(parent = environment(promise_ii_trial))
  controlled$promise_ii_predict <- function(...) c(current = 0, maximum = 0)
  controlled$promise_ii_trial <- promise_ii_trial
  environment(controlled$promise_ii_trial) <- controlled
  simulate <- promise_ii_simulate
  environment(simulate) <- controlled
  out <- simulate(.70, 2, 20, seed = 231, design = design)
  expect_true(all(out$trials$reason == "futile"))
  expect_true(all(out$trials$n == 60))
  expect_true(all(out$trials$success))
})

test_that("PROMISE II simulations reproduce across worker counts and preserve RNG state", {
  set.seed(99)
  before <- .Random.seed
  one <- promise_ii_simulate(.68, 4, 40, seed = 335, ncores = 1)
  expect_identical(.Random.seed, before)
  two <- promise_ii_simulate(.68, 4, 40, seed = 335, ncores = 2)
  expect_identical(one$trials, two$trials)
  expect_identical(.Random.seed, before)
})
