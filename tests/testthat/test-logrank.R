test_that("native log-rank checks lengths before accessing event indicators", {
  args <- list(
    groupa = c(1, 2), groupb = c(3, 4),
    groupacensored = c(1, 1), groupbcensored = c(1, 1), onlyz = FALSE
  )

  for (group in c("groupacensored", "groupbcensored")) {
    for (events in list(numeric(), 1, c(1, 1, 1))) {
      bad <- args
      bad[[group]] <- events
      for (onlyz in c(FALSE, TRUE)) {
        bad$onlyz <- onlyz
        expect_error(do.call(logrank_instance, bad), "same length in each group")
      }
    }
  }
})

test_that("native log-rank rejects invalid times before sorting", {
  args <- list(
    groupa = c(1, 2), groupb = c(3, 4),
    groupacensored = c(1, 1), groupbcensored = c(1, 1), onlyz = FALSE
  )

  for (group in c("groupa", "groupb")) {
    for (time in c(-1, NA_real_, NaN, Inf, -Inf)) {
      bad <- args
      bad[[group]][1] <- time
      expect_error(do.call(logrank_instance, bad), "times must be finite and non-negative")
    }
  }
})

test_that("native log-rank checks event values before integer conversion", {
  args <- list(
    groupa = c(1, 2), groupb = c(3, 4),
    groupacensored = c(1, 1), groupbcensored = c(1, 1), onlyz = FALSE
  )

  for (group in c("groupacensored", "groupbcensored")) {
    for (event in c(-1, 2, 0.5, 1.5, NA_real_, NaN, Inf, -Inf)) {
      bad <- args
      bad[[group]][1] <- event
      expect_error(do.call(logrank_instance, bad), "event indicators must be 0 or 1")
    }
  }
})

test_that("native log-rank requires one non-missing logical onlyz value", {
  for (onlyz in list(logical(), c(TRUE, FALSE), NA, 0, 1, NA_real_, "FALSE", NULL)) {
    expect_error(
      logrank_instance(c(1, 2), c(3, 4), c(1, 1), c(1, 1), onlyz),
      "'onlyz' must be a single non-missing logical value"
    )
  }
})

test_that("non-estimable log-rank comparisons do not return a valid p-value", {
  cases <- list(
    no_events = list(c(1, 2), c(3, 4), c(0, 0), c(0, 0)),
    all_events_tied = list(c(1, 1), c(1, 1), c(1, 1), c(1, 1)),
    no_overlap_at_events = list(c(1, 2), c(3, 4), c(0, 0), c(1, 1)),
    empty_a = list(numeric(), c(1, 2), numeric(), c(1, 1)),
    empty_b = list(c(1, 2), numeric(), c(1, 1), numeric()),
    both_empty = rep(list(numeric()), 4)
  )

  for (case in cases) {
    full <- do.call(logrank_test, case)
    expect_length(full, 3L)
    expect_true(all(is.nan(full)))
    expect_error(assert_logrank_estimable(full), "Log-rank analysis is non-estimable")
    z <- do.call(logrank_test, c(case, list(onlyz = TRUE)))
    expect_length(z, 1L)
    expect_true(is.nan(z))
  }
})

test_that("log-rank retains very small nonzero upper-tail probabilities", {
  for (n in c(40L, 120L)) {
    full <- logrank_test(seq_len(n), n + seq_len(n), rep(1L, n), rep(1L, n))
    expected <- pchisq(full[1], df = 1, lower.tail = FALSE)
    expect_gt(expected, 0)
    expect_lt(expected, .Machine$double.eps)
    expect_gt(full[3], 0)
    expect_equal(full[3] / expected, 1, tolerance = 1e-10)
  }
})

test_that("log-rank agrees with survdiff for ties, censoring, and unequal groups", {
  cases <- list(
    tied_and_unsorted = data.frame(
      time = c(6, 6, 5, 4, 3, 1, 6, 5, 5, 3, 2, 0),
      event = c(0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1),
      arm = rep(0:1, each = 6)
    ),
    no_events_in_one_arm = data.frame(
      time = c(10, 5, 4, 1, 3),
      event = c(0, 0, 0, 1, 1),
      arm = c(0, 0, 0, 1, 1)
    ),
    one_participant_per_arm = data.frame(time = c(1, 2), event = 1, arm = 0:1)
  )
  set.seed(42130)
  for (i in seq_len(25)) {
    n <- sample(4:35, 2, replace = TRUE)
    data <- data.frame(
      time = sample(0:20, sum(n), replace = TRUE),
      event = rbinom(sum(n), 1, 0.6),
      arm = rep(0:1, times = n)
    )
    # Ensure an event while both groups remain at risk.
    data$time[1] <- 0
    data$event[1] <- 1
    data$time[n[1] + 1] <- 21
    cases[[paste0("generated_", i)]] <- data
  }

  for (data in cases) {
    # Integer times give the same ties with or without near-tie adjustment.
    reference <- survival::survdiff(
      survival::Surv(time, event) ~ arm, data = data
    )
    a <- data[data$arm == 0, ]
    b <- data[data$arm == 1, ]
    full <- logrank_test(a$time, b$time, a$event, b$event)
    expected_z <- unname((reference$obs[1] - reference$exp[1]) / sqrt(reference$var[1, 1]))
    expect_equal(full[1], unname(reference$chisq), tolerance = 1e-10)
    expect_equal(full[2], expected_z, tolerance = 1e-10)
    expect_equal(full[3], pchisq(reference$chisq, 1, lower.tail = FALSE), tolerance = 1e-10)
    expect_equal(
      logrank_test(b$time, a$time, b$event, a$event),
      full * c(1, -1, 1), tolerance = 1e-10
    )
    expect_equal(
      logrank_test(rev(a$time), rev(b$time), rev(a$event), rev(b$event)),
      full, tolerance = 1e-10
    )
  }
})
