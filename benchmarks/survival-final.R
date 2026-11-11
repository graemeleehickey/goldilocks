# Bayesian survival final-analysis benchmark.
# Author: Graeme L. Hickey
# Run from the package root:
# Rscript benchmarks/survival-final.R 100 500 20 benchmarks/survival-final
# Arguments: legacy precision repetitions, direct precision repetitions,
# timing repetitions, and output prefix. No package defaults are changed.

pkgload::load_all(".", quiet = TRUE)
cli <- commandArgs(trailingOnly = TRUE)
legacy_repetitions <- if (length(cli) >= 1L) as.integer(cli[1]) else 100L
direct_repetitions <- if (length(cli) >= 2L) as.integer(cli[2]) else 500L
timing_repetitions <- if (length(cli) >= 3L) as.integer(cli[3]) else 20L
output <- if (length(cli) >= 4L) cli[4] else "benchmarks/survival-final"
stopifnot(
  legacy_repetitions >= 2L,
  direct_repetitions >= 2L,
  timing_repetitions >= 2L
)
reference_draws <- 20000000L
reference_batch <- 100000L
thresholds <- c(0.95, 0.975)
direct_draws <- c(1000L, 10000L, 50000L, 500000L)
legacy_commit <- "7cc2294"
legacy_source <- system2(
  "git",
  c("show", paste0(legacy_commit, ":R/analysis_final.R")),
  stdout = TRUE
)
if (!is.null(attr(legacy_source, "status"))) {
  stop("The pre-simplification commit must be available in this checkout.")
}
legacy_environment <- new.env(parent = asNamespace("goldilocks"))
eval(parse(text = legacy_source), legacy_environment)
legacy_final <- legacy_environment$analyse_final
direct_final <- goldilocks:::analyse_final

# A fixed dataset is used for each workload. Two constant-hazard datasets are
# deliberately selected near decision boundaries; this is not an OC study.
scenarios <- list(
  constant_10 = list(
    label = "Constant, 10% dropout",
    n = 400L,
    cutpoints = NULL,
    treatment = -log(0.78) / 12,
    control = -log(0.70) / 12,
    dropout = 0.10,
    h0 = 0,
    seed = 68100L,
    target = c(0.943, 0.947)
  ),
  constant_30 = list(
    label = "Constant, 30% dropout",
    n = 400L,
    cutpoints = NULL,
    treatment = -log(0.78) / 12,
    control = -log(0.70) / 12,
    dropout = 0.30,
    h0 = 0,
    seed = 68200L,
    target = c(0.978, 0.982)
  ),
  piecewise_30 = list(
    label = "Three intervals, 30% dropout",
    n = 400L,
    cutpoints = c(3, 6),
    treatment = c(0.015, 0.025, 0.035),
    control = c(0.025, 0.035, 0.045),
    dropout = 0.30,
    h0 = 0,
    seed = 681L
  ),
  unequal_dropout = list(
    label = "Three intervals, unequal dropout",
    n = 400L,
    cutpoints = c(3, 6),
    treatment = c(0.015, 0.025, 0.035),
    control = c(0.025, 0.035, 0.045),
    dropout = c(control = 0.10, treatment = 0.40),
    h0 = 0,
    seed = 68400L
  ),
  single_arm = list(
    label = "Single arm, 20% dropout",
    n = 120L,
    cutpoints = NULL,
    treatment = -log(0.70) / 12,
    control = NULL,
    dropout = 0.20,
    h0 = 0.40,
    seed = 68500L
  ),
  sparse = list(
    label = "Small sparse trial, 50% dropout",
    n = 80L,
    cutpoints = c(3, 6),
    treatment = c(0.012, 0.008, 0.006),
    control = c(0.018, 0.012, 0.009),
    dropout = 0.50,
    h0 = 0,
    seed = 68600L
  )
)

make_data <- function(scenario, seed) {
  set.seed(seed)
  data <- sim_trial_data(
    hazard_treatment = scenario$treatment,
    hazard_control = scenario$control,
    generation_cutpoints = scenario$cutpoints,
    N_total = scenario$n,
    lambda = 20,
    end_of_study = 12,
    block = 2,
    rand_ratio = c(control = 1, treatment = 1),
    prop_loss = scenario$dropout
  )
  data$subject_impute_success <- data$event == 0 & data$time < 12
  data
}

make_args <- function(scenario, data) {
  list(
    data_in = data,
    cutpoints = scenario$cutpoints,
    prior_surv_final = goldilocks:::normalize_gamma_prior(
      c(0.1, 0.1),
      length(scenario$cutpoints) + 1L,
      is.null(scenario$control),
      "prior_surv_final"
    ),
    N_mcmc = 1000L,
    single_arm = is.null(scenario$control),
    imputed_final = TRUE,
    method = "bayes-surv",
    N_impute = 500L,
    alternative = "less",
    h0 = scenario$h0,
    prior_bin = c(1, 1),
    bin_method = "mc",
    binary_imputation = "event-time",
    empty_interval = "prior",
    end_of_study = 12
  )
}

posterior_parameters <- function(args) {
  stats <- goldilocks:::post_sufficient_stats(
    args$data_in,
    args$cutpoints,
    args$single_arm
  )
  stats$shape <- 0.1 + stats$tot_events
  stats$rate <- 0.1 + stats$tot_time
  stats$width <- rep(
    diff(c(0, args$cutpoints, 12)),
    length(unique(stats$treatment))
  )
  stats
}

exact_probability <- function(args, stats) {
  if (!is.null(args$cutpoints)) {
    return(NA_real_)
  }
  treatment <- stats[stats$treatment == 1, ]
  if (args$single_arm) {
    return(pgamma(
      -log1p(-args$h0) / 12,
      shape = treatment$shape,
      rate = treatment$rate
    ))
  }
  control <- stats[stats$treatment == 0, ]
  pbeta(
    treatment$rate / (treatment$rate + control$rate),
    shape1 = treatment$shape,
    shape2 = control$shape
  )
}

effect_moments <- function(stats, single_arm) {
  survival_moments <- lapply(split(stats, stats$treatment), function(arm) {
    mean <- prod((arm$rate / (arm$rate + arm$width))^arm$shape)
    variance <- prod((arm$rate / (arm$rate + 2 * arm$width))^arm$shape) - mean^2
    c(mean = mean, variance = variance)
  })
  treatment <- survival_moments[["1"]]
  if (single_arm) {
    return(c(mean = 1 - unname(treatment[1]), variance = unname(treatment[2])))
  }
  control <- survival_moments[["0"]]
  c(
    mean = unname(control[1] - treatment[1]),
    variance = unname(control[2] + treatment[2])
  )
}

reference_probability <- function(args, stats, seed) {
  exact <- exact_probability(args, stats)
  if (is.finite(exact)) {
    return(list(
      probability = exact,
      mcse = 0,
      draws = 0L,
      method = "Gamma CDF / scaled-Gamma Beta identity"
    ))
  }
  # Independent batches bound memory. These draws are not reused in either
  # the repeated estimates or the runtime measurements.
  set.seed(seed)
  count <- 0
  for (batch in seq_len(reference_draws / reference_batch)) {
    cumulative <- lapply(split(stats, stats$treatment), function(arm) {
      total <- numeric(reference_batch)
      for (j in seq_len(nrow(arm))) {
        total <- total +
          arm$width[j] *
            rgamma(reference_batch, shape = arm$shape[j], rate = arm$rate[j])
      }
      total
    })
    count <- count + sum(cumulative[["1"]] < cumulative[["0"]])
  }
  probability <- count / reference_draws
  list(
    probability = probability,
    mcse = sqrt(probability * (1 - probability) / reference_draws),
    draws = reference_draws,
    method = "Independent Gamma draws"
  )
}

repeat_estimates <- function(calculator, args, repetitions, seed, label) {
  estimates <- matrix(
    NA_real_,
    repetitions,
    2L,
    dimnames = list(NULL, c("probability", "effect"))
  )
  seeds <- seed + seq_len(repetitions)
  for (i in seq_len(repetitions)) {
    set.seed(seeds[i])
    estimates[i, ] <- do.call(calculator, args)
    if (i %% 25L == 0L && label == "legacy") {
      cat("  legacy precision:", i, "/", repetitions, "\n")
      flush.console()
    }
  }
  data.frame(repetition = seq_len(repetitions), seed = seeds, estimates)
}

time_calculation <- function(calculator, args, seed) {
  set.seed(seed)
  invisible(do.call(calculator, args))
  gc()
  timing <- bench::mark(
    do.call(calculator, args),
    iterations = timing_repetitions,
    check = FALSE,
    memory = FALSE,
    filter_gc = FALSE
  )
  as.numeric(timing$time[[1]])
}

results <- list()
summary_rows <- list()
metadata <- list(
  author = "Graeme L. Hickey",
  date = "2026-11-11",
  source_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  legacy_commit = legacy_commit,
  legacy_source = legacy_source,
  current_final_source = deparse(direct_final),
  session = sessionInfo(),
  system = Sys.info(),
  legacy_repetitions = legacy_repetitions,
  direct_repetitions = direct_repetitions,
  timing_repetitions = timing_repetitions,
  reference_draws = reference_draws,
  scenarios = scenarios,
  thresholds = thresholds
)

for (index in seq_along(scenarios)) {
  name <- names(scenarios)[index]
  scenario <- scenarios[[index]]
  seed <- scenario$seed
  data <- make_data(scenario, seed)
  args <- make_args(scenario, data)
  if (!is.null(scenario$target)) {
    for (attempt in seq_len(1000L)) {
      seed <- scenario$seed + attempt - 1L
      data <- make_data(scenario, seed)
      args <- make_args(scenario, data)
      probability <- exact_probability(args, posterior_parameters(args))
      if (
        probability > scenario$target[1] && probability < scenario$target[2]
      ) {
        break
      }
    }
    stopifnot(
      probability > scenario$target[1],
      probability < scenario$target[2]
    )
  }
  stopifnot(any(data$subject_impute_success))
  stats <- posterior_parameters(args)
  reference <- reference_probability(args, stats, 700000L + index)
  moments <- effect_moments(stats, args$single_arm)
  cat(
    name,
    ": N =",
    nrow(data),
    "; dropout =",
    sum(data$loss_to_fu),
    "; reference probability =",
    format(reference$probability, digits = 8),
    "\n"
  )
  flush.console()

  # The current package's retained imputation path with 'error' is exactly the
  # former 'prior' calculation when every interval has exposure throughout.
  if (all(stats$n > 0)) {
    check_args <- args
    check_args$N_impute <- 3L
    set.seed(710000L + index)
    before <- do.call(legacy_final, check_args)
    before_rng <- .Random.seed
    check_args$empty_interval <- "error"
    set.seed(710000L + index)
    after <- do.call(direct_final, check_args)
    stopifnot(identical(before, after), identical(before_rng, .Random.seed))
  }

  legacy_estimates <- repeat_estimates(
    legacy_final,
    args,
    legacy_repetitions,
    1000000L + index * 10000L,
    "legacy"
  )
  direct_estimates <- repeat_estimates(
    direct_final,
    args,
    direct_repetitions,
    2000000L + index * 10000L,
    "direct"
  )
  legacy_sd <- sd(legacy_estimates$probability)
  matched_draws <- as.integer(ceiling(
    reference$probability *
      (1 - reference$probability) /
      legacy_sd^2
  ))
  methods <- c(
    list(legacy = c(impute = 500L, draws = 1000L)),
    setNames(
      lapply(direct_draws, function(n) c(impute = 0L, draws = n)),
      paste0("direct_", direct_draws)
    ),
    list(direct_matched = c(impute = 0L, draws = matched_draws))
  )
  timings <- list()
  for (method in names(methods)) {
    settings <- methods[[method]]
    timing_args <- args
    timing_args$N_mcmc <- settings["draws"]
    calculator <- if (method == "legacy") legacy_final else direct_final
    timings[[method]] <- time_calculation(
      calculator,
      timing_args,
      3000000L + index * 10000L + match(method, names(methods))
    )
    estimates <- if (method == "legacy") {
      legacy_estimates
    } else {
      if (method == "direct_1000") direct_estimates else NULL
    }
    empirical_sd <- if (!is.null(estimates)) {
      sd(estimates$probability)
    } else {
      NA_real_
    }
    theoretical_sd <- if (method != "legacy") {
      sqrt(
        reference$probability * (1 - reference$probability) / settings["draws"]
      )
    } else {
      NA_real_
    }
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      scenario = name,
      label = scenario$label,
      method = method,
      n = scenario$n,
      intervals = length(scenario$cutpoints) + 1L,
      single_arm = args$single_arm,
      data_seed = seed,
      observed_events = sum(data$event),
      observed_dropouts = sum(data$loss_to_fu),
      reference_probability = reference$probability,
      reference_mcse = reference$mcse,
      reference_draws = reference$draws,
      reference_effect = moments["mean"],
      n_impute = settings["impute"],
      n_mcmc = settings["draws"],
      median_ms = 1000 * median(timings[[method]]),
      q25_ms = 1000 * unname(quantile(timings[[method]], 0.25)),
      q75_ms = 1000 * unname(quantile(timings[[method]], 0.75)),
      empirical_repetitions = if (is.null(estimates)) 0L else nrow(estimates),
      empirical_probability_sd = empirical_sd,
      theoretical_probability_sd = theoretical_sd,
      probability_bias = if (!is.null(estimates)) {
        mean(estimates$probability) - reference$probability
      } else {
        NA_real_
      },
      empirical_effect_sd = if (!is.null(estimates)) {
        sd(estimates$effect)
      } else {
        NA_real_
      },
      theoretical_effect_sd = if (method != "legacy") {
        sqrt(moments["variance"] / settings["draws"])
      } else {
        NA_real_
      },
      effect_bias = if (!is.null(estimates)) {
        mean(estimates$effect) - moments["mean"]
      } else {
        NA_real_
      },
      sd_lower = if (!is.null(estimates)) {
        empirical_sd *
          sqrt(
            (nrow(estimates) - 1) /
              qchisq(0.975, nrow(estimates) - 1)
          )
      } else {
        NA_real_
      },
      sd_upper = if (!is.null(estimates)) {
        empirical_sd *
          sqrt(
            (nrow(estimates) - 1) /
              qchisq(0.025, nrow(estimates) - 1)
          )
      } else {
        NA_real_
      },
      row.names = NULL
    )
  }
  results[[name]] <- list(
    scenario = scenario,
    data_seed = seed,
    data = data,
    args = args,
    sufficient_statistics = stats,
    reference = reference,
    effect_moments = moments,
    legacy_estimates = legacy_estimates,
    direct_estimates = direct_estimates,
    timings_seconds = timings
  )
  summary <- do.call(rbind, summary_rows)
  saveRDS(
    list(metadata = metadata, results = results, summary = summary),
    paste0(output, "-results.rds")
  )
  write.csv(summary, paste0(output, "-results.csv"), row.names = FALSE)
  cat(
    "  completed; legacy SD =",
    format(legacy_sd, digits = 4),
    "; equivalent direct draws =",
    matched_draws,
    "\n"
  )
  flush.console()
}

# Exact binomial crossing probabilities for the direct estimator, including
# the package's strict success comparison. These describe numerical decision
# error conditional on a fixed posterior probability, not trial type I error.
decisions <- expand.grid(
  threshold = thresholds,
  draws = direct_draws,
  gap = c(-0.01, -0.005, -0.002, -0.001, 0.001, 0.002, 0.005, 0.01)
)
decisions$probability <- decisions$threshold + decisions$gap
decisions$success_probability <- pbinom(
  floor(decisions$threshold * decisions$draws),
  size = decisions$draws,
  prob = decisions$probability,
  lower.tail = FALSE
)
decisions$wrong_decision_probability <- ifelse(
  decisions$gap > 0,
  1 - decisions$success_probability,
  decisions$success_probability
)
write.csv(decisions, paste0(output, "-decisions.csv"), row.names = FALSE)
saveRDS(
  list(
    metadata = metadata,
    results = results,
    summary = summary,
    decisions = decisions
  ),
  paste0(output, "-results.rds")
)
cat("Saved benchmark:", output, "\n")
