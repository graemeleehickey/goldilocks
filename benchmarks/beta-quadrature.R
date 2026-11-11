# Bayesian binary quadrature benchmark.
# Author: Graeme L. Hickey
# Run from the package root:
# Rscript benchmarks/beta-quadrature.R 200 5 benchmarks/beta-quadrature-results.csv
# Arguments: probability iterations, trial-batch iterations, output CSV.
# Experimental calculations replace the loaded R function only during timing;
# no package source files are changed.

for (package in c("bench", "pkgload", "Rmpfr")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Install '", package, "' to run this benchmark.")
  }
}
pkgload::load_all(".", quiet = TRUE, compile = FALSE)
arguments <- commandArgs(trailingOnly = TRUE)
iterations <- if (length(arguments)) as.integer(arguments[1]) else 200L
trial_iterations <- if (length(arguments) > 1) as.integer(arguments[2]) else 5L
output <- if (length(arguments) > 2) arguments[3] else {
  "benchmarks/beta-quadrature-results.csv"
}
stopifnot(is.finite(iterations), iterations > 0,
          is.finite(trial_iterations), trial_iterations > 0)

original <- getFromNamespace("analyse_bayes_bin_from_counts", "goldilocks")

quantile_probability <- function(treatment, control, h0, alternative,
                                 split = FALSE) {
  integrand <- function(u) {
    threshold <- stats::qbeta(u, treatment$alpha, treatment$beta) - h0
    stats::pbeta(threshold, control$alpha, control$beta,
                 lower.tail = alternative == "greater")
  }
  breaks <- c(0, 1)
  if (split) {
    # Locate the change in the control CDF on the treatment quantile scale.
    control_quantiles <- stats::qbeta(
      c(1e-10, 0.01, 0.5, 0.99, 1 - 1e-10),
      control$alpha, control$beta
    )
    breaks <- sort(unique(c(breaks, stats::pbeta(
      control_quantiles + h0, treatment$alpha, treatment$beta
    ))))
  }
  intervals <- length(breaks) - 1L
  probability <- sum(vapply(seq_len(intervals), function(i) {
    width <- breaks[i + 1L] - breaks[i]
    # The integrand is bounded by [0, 1]. On a negligible interval, midpoint
    # evaluation has absolute error at most its width and avoids roundoff in
    # adaptive integration extremely close to u = 0 or u = 1.
    if (width <= 1e-10 / intervals) {
      return(width * integrand(breaks[i] + width / 2))
    }
    stats::integrate(
      integrand, breaks[i], breaks[i + 1L], subdivisions = 200L,
      rel.tol = 1e-10, abs.tol = 1e-10 / intervals
    )$value
  }, numeric(1)))
  min(1, max(0, probability))
}

make_calculator <- function(split) {
  calculator <- original
  environment(calculator) <- new.env(parent = environment(original))
  environment(calculator)$benchmark_probability <- function(treatment, control,
                                                           h0, alternative) {
    quantile_probability(treatment, control, h0, alternative, split = split)
  }
  expressions <- as.list(body(calculator))
  location <- which(vapply(expressions, function(x) {
    is.call(x) && identical(x[[1]], as.name("<-")) &&
      identical(x[[2]], as.name("beta_binomial_difference_success"))
  }, logical(1)))
  stopifnot(length(location) == 1L)
  expressions[[location]] <- quote(
    beta_binomial_difference_success <- function(treatment, control) {
      benchmark_probability(treatment, control, h0, alternative)
    }
  )
  body(calculator) <- as.call(expressions)
  calculator
}

calculators <- list(
  current = original,
  quantile = make_calculator(FALSE),
  split_quantile = make_calculator(TRUE)
)

cases <- data.frame(
  case = c("small_alternative", "typical_null", "typical_alternative",
           "near_success_threshold", "nonzero_margin", "sparse_small",
           "jeffreys_equal", "jeffreys_sparse", "concentrated_zero",
           "concentrated_zero_large", "concentrated_all_events",
           "concentrated_rare_difference", "concentrated_central",
           "concentrated_central_alternative", "small_tail", "sharp_transition"),
  category = c(rep("typical", 5), rep("sparse", 3), rep("stress", 8)),
  events_treatment = c(10, 75, 55, 60, 55, 0, 0, 0, 0, 0, 10000, 0, 3000, 2900, 10, 0),
  n_treatment = c(50, 250, 250, 250, 250, 10, 10, 20, 10000, 100000, 10000, 10000, 10000, 10000, 250, 20),
  events_control = c(15, 75, 75, 75, 75, 1, 0, 1, 0, 0, 10000, 1, 3000, 3000, 75, 3000),
  n_control = c(50, 250, 250, 250, 250, 10, 10, 20, 10000, 100000, 10000, 10000, 10000, 10000, 250, 10000),
  prior_a = c(rep(1, 6), 0.5, 0.5, rep(1, 8)),
  prior_b = c(rep(1, 6), 0.5, 0.5, rep(1, 8)),
  h0 = c(rep(0, 4), -0.05, rep(0, 11))
)

case_arguments <- function(case, alternative) {
  list(events_control = case$events_control, n_control = case$n_control,
       events_treatment = case$events_treatment, n_treatment = case$n_treatment,
       single_arm = FALSE, alternative = alternative, h0 = case$h0,
       prior_bin = c(case$prior_a, case$prior_b),
       bin_method = "quadrature", N_mcmc = NULL)
}

# For zero margins and integer Beta shapes, integrate the finite binomial
# expansion of the Beta CDF using 128-bit arithmetic. The positive-term sum
# gives P(T > C) independently of numerical quadrature.
integer_reference <- function(at, bt, ac, bc, alternative) {
  i <- Rmpfr::mpfr(seq.int(0, at - 1), 128)
  bt <- Rmpfr::mpfr(bt, 128)
  ac <- Rmpfr::mpfr(ac, 128)
  bc <- Rmpfr::mpfr(bc, 128)
  greater <- sum(exp(
    Rmpfr::lbeta(ac + i, bt + bc) - log(bt + i) -
      Rmpfr::lbeta(1 + i, bt) - Rmpfr::lbeta(ac, bc)
  ))
  as.numeric(if (alternative == "greater") greater else 1 - greater)
}

# For noninteger shapes or nonzero margins, use tightly partitioned density
# integration with x = sin(theta)^2 to remove Jeffreys-prior endpoint
# singularities. Verify agreement after swapping which arm is integrated.
density_reference <- function(at, bt, ac, bc, h0, alternative) {
  probabilities <- c(0, 1e-12, 1e-8, 1e-4, 0.01, 0.1, 0.5,
                     0.9, 0.99, 1 - 1e-4, 1 - 1e-8, 1 - 1e-12, 1)
  breaks <- sort(unique(c(
    stats::qbeta(probabilities, at, bt),
    pmin(1, pmax(0, stats::qbeta(probabilities, ac, bc) + h0))
  )))
  breaks <- sort(unique(asin(sqrt(breaks))))
  sum(vapply(seq_len(length(breaks) - 1L), function(i) {
    stats::integrate(function(theta) {
      x <- sin(theta)^2
      2 * sin(theta) * cos(theta) * stats::dbeta(x, at, bt) * stats::pbeta(
        x - h0, ac, bc, lower.tail = alternative == "greater"
      )
    }, breaks[i], breaks[i + 1L], rel.tol = 1e-12, abs.tol = 1e-13,
    subdivisions = 1000L)$value
  }, numeric(1)))
}

reference_probability <- function(case, alternative) {
  at <- case$events_treatment + case$prior_a
  bt <- case$n_treatment - case$events_treatment + case$prior_b
  ac <- case$events_control + case$prior_a
  bc <- case$n_control - case$events_control + case$prior_b
  if (case$h0 == 0 && at == ac && bt == bc) {
    return(list(value = 0.5, method = "symmetry"))
  }
  if (case$h0 == 0 && all(c(at, bt, ac, bc) %% 1 == 0)) {
    return(list(value = integer_reference(at, bt, ac, bc, alternative),
                method = "128-bit finite sum"))
  }
  value <- density_reference(at, bt, ac, bc, case$h0, alternative)
  other <- density_reference(ac, bc, at, bt, -case$h0,
                             if (alternative == "less") "greater" else "less")
  stopifnot(abs(value - other) < 1e-10)
  list(value = value, method = "partitioned density; both arm orientations")
}

run_time <- format(Sys.time(), tz = "UTC", usetz = TRUE)
source_commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
metadata <- list(run_utc = run_time, source_commit = source_commit,
                 R_version = as.character(getRversion()),
                 package_version = as.character(utils::packageVersion("goldilocks")),
                 bench_version = as.character(utils::packageVersion("bench")),
                 Rmpfr_version = as.character(utils::packageVersion("Rmpfr")),
                 platform = R.version$platform)
rows <- list()
timing_columns <- function(timing) {
  list(median_seconds = as.numeric(timing$median),
       minimum_seconds = as.numeric(timing$min),
       q25_seconds = as.numeric(stats::quantile(timing$time[[1]], 0.25)),
       q75_seconds = as.numeric(stats::quantile(timing$time[[1]], 0.75)),
       n_iterations = timing$n_itr, n_gc = timing$n_gc)
}

set.seed(92501)
for (index in seq_len(nrow(cases))) {
  case <- cases[index, ]
  message("Probability case ", index, "/", nrow(cases), ": ", case$case)
  for (alternative in c("less", "greater")) {
    input <- case_arguments(case, alternative)
    reference <- reference_probability(case, alternative)
    # Vary ordering to reduce systematic timing advantages from running first.
    for (variant in sample(names(calculators))) {
      calculator <- calculators[[variant]]
      value <- do.call(calculator, input)$success
      timing <- bench::mark(
        do.call(calculator, input), iterations = iterations,
        check = FALSE, memory = FALSE, filter_gc = TRUE
      )
      rows[[length(rows) + 1L]] <- as.data.frame(c(
        list(workload = "probability"), as.list(case),
        list(alternative = alternative, variant = variant, probability = value,
             reference = reference$value, absolute_error = abs(value - reference$value),
             reference_method = reference$method), timing_columns(timing), metadata
      ))
    }
  }
}

probability_results <- do.call(rbind, rows)
split_errors <- subset(probability_results, variant == "split_quantile")$absolute_error
stopifnot(all(split_errors < 1e-8))

with_calculator <- function(calculator, code) {
  namespace <- asNamespace("goldilocks")
  name <- "analyse_bayes_bin_from_counts"
  saved <- get(name, envir = namespace)
  on.exit({
    unlockBinding(name, namespace)
    assign(name, saved, envir = namespace)
    lockBinding(name, namespace)
  })
  unlockBinding(name, namespace)
  assign(name, calculator, envir = namespace)
  lockBinding(name, namespace)
  force(code)
}

simulate_batch <- function() {
  goldilocks::sim_trials(
    hazard_control = -log(1 - 0.30) / 12,
    hazard_treatment = -log(1 - 0.22) / 12,
    N_total = 500, interim_look = c(100, 200, 300, 400),
    end_of_study = 12, lambda = 10, prior_surv = c(0.5, 20),
    prior_surv_final = c(0.5, 20), prior_bin = c(1, 1),
    alternative = "less", method = "bayes-bin", bin_method = "quadrature",
    N_impute = 100, N_mcmc = NULL, Fn = 0.05, Sn = 0.9, Qn = 1,
    prob_ha = 0.95, N_trials = 10, backend = "sequential", seed = 92502
  )
}
baseline <- simulate_batch()
decision_columns <- c("N_enrolled", "stopping_reason", "trial_success")
for (variant in sample(names(calculators))) {
  message("Trial batch: ", variant)
  calculator <- calculators[[variant]]
  result <- suppressWarnings(with_calculator(calculator, simulate_batch()))
  failures <- nrow(result$failures)
  decisions_identical <- if (failures) NA else {
    identical(result$sims[decision_columns], baseline$sims[decision_columns])
  }
  if (variant == "split_quantile") stopifnot(failures == 0L, decisions_identical)
  if (failures) {
    # A run that omits trials is not a comparable full-workload timing.
    timing_result <- list(median_seconds = NA_real_, minimum_seconds = NA_real_,
                          q25_seconds = NA_real_, q75_seconds = NA_real_,
                          n_iterations = 0L, n_gc = NA_integer_)
    message(failures, "/10 trials failed; no full-workload timing reported")
  } else {
    timing <- with_calculator(calculator, bench::mark(
      simulate_batch(), iterations = trial_iterations,
      check = FALSE, memory = FALSE, filter_gc = FALSE
    ))
    timing_result <- timing_columns(timing)
  }
  rows[[length(rows) + 1L]] <- as.data.frame(c(
    list(workload = "10 adaptive trials", case = "N500_four_interims_100_imputations",
         category = "trial", alternative = "less", variant = variant,
         probability = if (failures) NA_real_ else mean(result$sims$trial_success),
         reference = mean(baseline$sims$trial_success),
         absolute_error = if (failures) NA_real_ else {
           max(abs(result$sims$post_prob_ha - baseline$sims$post_prob_ha))
         },
         reference_method = "matched simulated histories; compare final scores",
         n_requested = 10L, n_analyzed = nrow(result$sims), n_failed = failures,
         decisions_identical = decisions_identical,
         error_message = if (failures) result$failures$message[1] else ""),
    timing_result, metadata
  ))
}

results <- dplyr::bind_rows(rows)
results <- results |>
  dplyr::group_by(workload, case, alternative) |>
  dplyr::mutate(time_ratio_to_current = median_seconds / median_seconds[variant == "current"]) |>
  dplyr::ungroup()
utils::write.csv(results, output, row.names = FALSE, na = "")
print(results |>
  dplyr::group_by(workload, category, variant) |>
  dplyr::summarise(median_seconds = median(median_seconds),
                   median_time_ratio = median(time_ratio_to_current),
                   maximum_absolute_error = if (all(is.na(absolute_error))) {
                     NA_real_
                   } else max(absolute_error, na.rm = TRUE), .groups = "drop"),
  n = Inf)
cat("\nResults saved to", output, "\n")
