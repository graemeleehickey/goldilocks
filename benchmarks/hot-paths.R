# Optional maintainer benchmarks for simulation hot paths.
#
# Run from the package root with:
#   source("benchmarks/hot-paths.R")

if (!requireNamespace("bench", quietly = TRUE)) {
  stop("Install the 'bench' package to run these benchmarks.", call. = FALSE)
}

load_goldilocks <- function() {
  if (
    requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")
  ) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    library(goldilocks)
  }
}

load_goldilocks()

ppwe_internal <- getFromNamespace("ppwe", "goldilocks")
haz_to_prop_internal <- getFromNamespace("haz_to_prop", "goldilocks")
bayes_surv_effect_draws_kernel_internal <- getFromNamespace(
  "bayes_surv_effect_draws_kernel",
  "goldilocks"
)
endpoint_interval_widths_internal <- getFromNamespace(
  "endpoint_interval_widths",
  "goldilocks"
)
posterior_internal <- getFromNamespace("posterior", "goldilocks")
posterior_sufficient_stats_internal <- getFromNamespace(
  "posterior_sufficient_stats",
  "goldilocks"
)
posterior_from_sufficient_stats_internal <- getFromNamespace(
  "posterior_from_sufficient_stats",
  "goldilocks"
)
posterior_from_sufficient_stats_kernel_internal <- getFromNamespace(
  "posterior_from_sufficient_stats_kernel",
  "goldilocks"
)
normalize_gamma_prior_internal <- getFromNamespace(
  "normalize_gamma_prior",
  "goldilocks"
)
impute_data_internal <- getFromNamespace("impute_data", "goldilocks")
impute_predictive_draws_internal <- getFromNamespace(
  "impute_predictive_draws",
  "goldilocks"
)
materialize_predictive_draw_internal <- getFromNamespace(
  "materialize_predictive_draw",
  "goldilocks"
)
analyse_data_internal <- getFromNamespace("analyse_data", "goldilocks")
analyse_bayes_surv_sufficient_stats_kernel_internal <- getFromNamespace(
  "analyse_bayes_surv_sufficient_stats_kernel",
  "goldilocks"
)
prepare_predictive_survival_state_internal <- getFromNamespace(
  "prepare_predictive_survival_state",
  "goldilocks"
)
test_stop_success_internal <- getFromNamespace(
  "test_stop_success",
  "goldilocks"
)
predictive_binary_count_states_internal <- getFromNamespace(
  "predictive_binary_count_states",
  "goldilocks"
)
analyse_predictive_binary_counts_internal <- getFromNamespace(
  "analyse_predictive_binary_counts",
  "goldilocks"
)
cumulative_hazard_to_probability_internal <- getFromNamespace(
  "cumulative_hazard_to_probability",
  "goldilocks"
)
probability_to_cumulative_hazard_internal <- getFromNamespace(
  "probability_to_cumulative_hazard",
  "goldilocks"
)
cox_wald_test_internal <- getFromNamespace("cox_wald_test", "goldilocks")
coxph_fit_compatibility_internal <- getFromNamespace(
  "coxph_fit_compatibility",
  "goldilocks"
)
cox_compatibility <- coxph_fit_compatibility_internal()
message(
  "Cox fast path: ",
  if (cox_compatibility$compatible) "enabled" else "disabled",
  " (",
  cox_compatibility$reason,
  ")"
)

iterations <- as.integer(Sys.getenv("GOLDILOCKS_BENCHMARK_ITERATIONS", "5"))
if (is.na(iterations) || iterations < 1L) {
  stop(
    "GOLDILOCKS_BENCHMARK_ITERATIONS must be a positive integer.",
    call. = FALSE
  )
}

set.seed(4242)

cutpoints_piecewise <- c(6, 12)
end_of_study <- 36
n_intervals_piecewise <- length(cutpoints_piecewise) + 1L
analysis_interval_widths <- endpoint_interval_widths_internal(
  cutpoints_piecewise,
  end_of_study
)

hazard_matrix <- matrix(
  stats::rgamma(15000, shape = 2, rate = 80),
  ncol = n_intervals_piecewise
)

cumulative_hazard_values <- c(
  0,
  10^seq(-20, 3, length.out = 100000),
  Inf
)
probability_values <- c(
  0,
  10^seq(-20, -1, length.out = 100000),
  1
)

posterior_draws <- array(
  stats::rgamma(30000, shape = 2, rate = 80),
  dim = c(5000, n_intervals_piecewise, 2)
)

posterior_data <- data.frame(
  time = pmin(stats::rexp(600, rate = 0.03), end_of_study),
  event = stats::rbinom(600, size = 1, prob = 0.65),
  treatment = rep(0:1, each = 300)
)
posterior_stats <- posterior_sufficient_stats_internal(
  data = posterior_data,
  cutpoints = cutpoints_piecewise,
  single_arm = FALSE
)
posterior_prior <- normalize_gamma_prior_internal(
  c(0.1, 0.1),
  n_intervals = n_intervals_piecewise,
  single_arm = FALSE,
  name = "prior_surv"
)

cox_data <- data.frame(
  time = round(stats::rexp(1000, rate = 0.03), digits = 1),
  event = stats::rbinom(1000, size = 1, prob = 0.65),
  treatment = rep(0:1, each = 500)
)

imputation_data <- data.frame(
  time = pmin(stats::rexp(500, rate = 0.025), end_of_study),
  treatment = rep(0:1, length.out = 500),
  event = stats::rbinom(500, size = 1, prob = 0.55),
  id = seq_len(500),
  subject_enrolled = c(rep(TRUE, 360), rep(FALSE, 140)),
  subject_impute_success = FALSE,
  subject_impute_futility = c(rep(FALSE, 360), rep(TRUE, 140))
)
imputation_data$subject_impute_success <-
  imputation_data$subject_enrolled &
  imputation_data$event == 0 &
  imputation_data$time < end_of_study

imputation_hazard <- array(
  c(0.018, 0.026, 0.034, 0.026, 0.034, 0.042),
  dim = c(1, n_intervals_piecewise, 2)
)
imputation_hazards <- array(
  stats::rgamma(100 * n_intervals_piecewise * 2, shape = 2, rate = 80),
  dim = c(100, n_intervals_piecewise, 2)
)

scalar_predictive_imputations <- function() {
  current <- vector("list", dim(imputation_hazards)[1])
  maximum <- vector("list", dim(imputation_hazards)[1])
  for (draw in seq_len(dim(imputation_hazards)[1])) {
    current[[draw]] <- impute_data_internal(
      data_in = imputation_data,
      hazard = imputation_hazards[draw, , , drop = FALSE],
      end_of_study = end_of_study,
      cutpoints = cutpoints_piecewise,
      type = "success",
      single_arm = FALSE
    )
    maximum[[draw]] <- impute_data_internal(
      data_in = current[[draw]],
      hazard = imputation_hazards[draw, , , drop = FALSE],
      end_of_study = end_of_study,
      cutpoints = cutpoints_piecewise,
      type = "futility",
      single_arm = FALSE
    )
  }
  list(current = current, maximum = maximum)
}

set.seed(1007)
binary_predictive_imputations <- impute_predictive_draws_internal(
  data_in = imputation_data,
  hazards = imputation_hazards,
  end_of_study = end_of_study,
  cutpoints = cutpoints_piecewise,
  single_arm = FALSE,
  binary_imputation = "bernoulli",
  check_futility = TRUE
)
binary_count_states <- predictive_binary_count_states_internal(
  data_in = imputation_data,
  imputations = binary_predictive_imputations,
  single_arm = FALSE,
  check_futility = TRUE
)

survival_analysis_state <- prepare_predictive_survival_state_internal(
  data_in = imputation_data,
  imputations = binary_predictive_imputations,
  check_futility = TRUE
)

materialized_survival_analyses <- function(method, N_mcmc) {
  output <- numeric(binary_predictive_imputations$n_draws)
  for (draw in seq_len(binary_predictive_imputations$n_draws)) {
    current <- materialize_predictive_draw_internal(
      data_in = imputation_data,
      imputations = binary_predictive_imputations,
      draw = draw
    )
    maximum <- materialize_predictive_draw_internal(
      data_in = current,
      imputations = binary_predictive_imputations,
      draw = draw,
      include_future = TRUE
    )

    if (method == "bayes-surv") {
      current_stats <- posterior_sufficient_stats_internal(
        data = current,
        cutpoints = cutpoints_piecewise,
        single_arm = FALSE,
        rows = current$subject_enrolled
      )
      maximum_stats <- posterior_sufficient_stats_internal(
        data = maximum,
        cutpoints = cutpoints_piecewise,
        single_arm = FALSE
      )
      current_result <- analyse_bayes_surv_sufficient_stats_kernel_internal(
        data_summ = current_stats,
        cutpoints = cutpoints_piecewise,
        end_of_study = end_of_study,
        prior_surv = posterior_prior,
        N_mcmc = N_mcmc,
        single_arm = FALSE,
        alternative = "less",
        h0 = 0,
        empty_interval = "prior",
        interval_widths = analysis_interval_widths
      )
      maximum_result <- analyse_bayes_surv_sufficient_stats_kernel_internal(
        data_summ = maximum_stats,
        cutpoints = cutpoints_piecewise,
        end_of_study = end_of_study,
        prior_surv = posterior_prior,
        N_mcmc = N_mcmc,
        single_arm = FALSE,
        alternative = "less",
        h0 = 0,
        empty_interval = "prior",
        interval_widths = analysis_interval_widths
      )
    } else {
      current_data <- current[
        current$subject_enrolled,
        c("time", "event", "treatment"),
        drop = FALSE
      ]
      maximum_data <- maximum[,
        c("time", "event", "treatment"),
        drop = FALSE
      ]
      current_result <- analyse_data_internal(
        data = current_data,
        cutpoints = cutpoints_piecewise,
        end_of_study = end_of_study,
        prior_surv = c(0.1, 0.1),
        N_mcmc = N_mcmc,
        single_arm = FALSE,
        method = method,
        alternative = "less",
        h0 = 0
      )
      maximum_result <- analyse_data_internal(
        data = maximum_data,
        cutpoints = cutpoints_piecewise,
        end_of_study = end_of_study,
        prior_surv = c(0.1, 0.1),
        N_mcmc = N_mcmc,
        single_arm = FALSE,
        method = method,
        alternative = "less",
        h0 = 0
      )
    }

    output[draw] <- current_result$success + maximum_result$success
  }
  output
}

vector_survival_analyses <- function(method, N_mcmc) {
  output <- numeric(binary_predictive_imputations$n_draws)
  for (draw in seq_len(binary_predictive_imputations$n_draws)) {
    result <- test_stop_success_internal(
      analysis_state = survival_analysis_state,
      imputations = binary_predictive_imputations,
      draw = draw,
      end_of_study = end_of_study,
      cutpoints = cutpoints_piecewise,
      interval_widths = if (method == "bayes-surv") {
        analysis_interval_widths
      } else {
        NULL
      },
      single_arm = FALSE,
      prior_surv = posterior_prior,
      N_mcmc = N_mcmc,
      method = method,
      alternative = "less",
      h0 = 0,
      empty_interval = "prior",
      check_futility = TRUE,
      prob_ha = 0.95,
      mc_conf_level = 0.95
    )
    output[draw] <- result$success_now$success + result$success_max$success
  }
  output
}

for (method in c("logrank", "cox")) {
  if (
    !identical(
      materialized_survival_analyses(method, 1L),
      vector_survival_analyses(method, 1L)
    )
  ) {
    stop(method, " vector benchmark does not match its materialized reference.")
  }
}
set.seed(1009)
materialized_bayes_reference <- materialized_survival_analyses(
  "bayes-surv",
  300L
)
set.seed(1009)
vector_bayes_reference <- vector_survival_analyses("bayes-surv", 300L)
if (!identical(materialized_bayes_reference, vector_bayes_reference)) {
  stop("Bayesian vector benchmark does not match its materialized reference.")
}

materialized_binary_analyses <- function() {
  current_successes <- 0L
  maximum_successes <- 0L
  for (draw in seq_len(binary_predictive_imputations$n_draws)) {
    current <- materialize_predictive_draw_internal(
      data_in = imputation_data,
      imputations = binary_predictive_imputations,
      draw = draw
    )
    current <- current[
      current$subject_enrolled,
      c("time", "event", "treatment"),
      drop = FALSE
    ]
    maximum <- materialize_predictive_draw_internal(
      data_in = imputation_data,
      imputations = binary_predictive_imputations,
      draw = draw,
      include_future = TRUE
    )
    maximum <- maximum[, c("time", "event", "treatment"), drop = FALSE]
    current_result <- analyse_data_internal(
      data = current,
      end_of_study = end_of_study,
      cutpoints = cutpoints_piecewise,
      single_arm = FALSE,
      prior_surv = c(0.1, 0.1),
      N_mcmc = 1L,
      method = "bayes-bin",
      alternative = "less",
      h0 = 0,
      prior_bin = c(1, 1),
      bin_method = "quadrature",
      empty_interval = "prior"
    )
    maximum_result <- analyse_data_internal(
      data = maximum,
      end_of_study = end_of_study,
      cutpoints = cutpoints_piecewise,
      single_arm = FALSE,
      prior_surv = c(0.1, 0.1),
      N_mcmc = 1L,
      method = "bayes-bin",
      alternative = "less",
      h0 = 0,
      prior_bin = c(1, 1),
      bin_method = "quadrature",
      empty_interval = "prior"
    )
    current_successes <- current_successes +
      as.integer(current_result$success > 0.95)
    maximum_successes <- maximum_successes +
      as.integer(maximum_result$success > 0.95)
  }
  c(current = current_successes, maximum = maximum_successes)
}

count_based_binary_analyses <- function() {
  result <- analyse_predictive_binary_counts_internal(
    count_states = binary_count_states,
    single_arm = FALSE,
    N_mcmc = 1L,
    method = "bayes-bin",
    alternative = "less",
    h0 = 0,
    prior_bin = c(1, 1),
    bin_method = "quadrature",
    check_futility = TRUE,
    prob_ha = 0.95,
    mc_conf_level = 0.95
  )
  c(
    current = result$current_successes,
    maximum = result$maximum_successes
  )
}

if (!identical(materialized_binary_analyses(), count_based_binary_analyses())) {
  stop("Binary count benchmark does not match its materialized reference.")
}
binary_reuse <- analyse_predictive_binary_counts_internal(
  count_states = binary_count_states,
  single_arm = FALSE,
  N_mcmc = 1L,
  method = "bayes-bin",
  alternative = "less",
  h0 = 0,
  prior_bin = c(1, 1),
  bin_method = "quadrature",
  check_futility = TRUE,
  prob_ha = 0.95,
  mc_conf_level = 0.95
)$reuse
message(
  "Binary count-state reuse: ",
  binary_reuse$cache_hits,
  "/",
  binary_reuse$analysis_requests,
  " analyses (",
  formatC(100 * binary_reuse$cache_hit_rate, digits = 1, format = "f"),
  "%)"
)

trial_cutpoints <- c(4, 8)
trial_hazard_control <- prop_to_haz(
  c(0.20, 0.30, 0.35),
  trial_cutpoints,
  end_of_study
)
trial_hazard_treatment <- prop_to_haz(
  c(0.10, 0.18, 0.22),
  trial_cutpoints,
  end_of_study
)

benchmark_results <- bench::mark(
  cumulative_hazard_to_probability = {
    cumulative_hazard_to_probability_internal(cumulative_hazard_values)
  },
  probability_to_cumulative_hazard = {
    probability_to_cumulative_hazard_internal(probability_values)
  },
  enrollment_constant = {
    set.seed(1000)
    enrollment(lambda = 0.001, N_total = 10000)
  },
  enrollment_piecewise = {
    set.seed(1000)
    enrollment(
      lambda = c(0.001, 0.05, 1),
      lambda_time = c(100, 250.5),
      N_total = 10000
    )
  },
  ppwe_piecewise = {
    ppwe_internal(
      hazard = hazard_matrix,
      end_of_study = end_of_study,
      cutpoints = cutpoints_piecewise
    )
  },
  haz_to_prop_piecewise = {
    haz_to_prop_internal(
      post = posterior_draws,
      cutpoints = cutpoints_piecewise,
      end_of_study = end_of_study,
      single_arm = FALSE
    )
  },
  bayes_surv_effect_kernel_piecewise = {
    bayes_surv_effect_draws_kernel_internal(
      post_lambda = posterior_draws,
      interval_widths = analysis_interval_widths,
      single_arm = FALSE
    )
  },
  posterior_piecewise = {
    set.seed(1001)
    posterior_internal(
      data = posterior_data,
      cutpoints = cutpoints_piecewise,
      prior_surv = c(0.1, 0.1),
      N_mcmc = 1000,
      single_arm = FALSE
    )
  },
  posterior_from_stats_checked = {
    set.seed(1001)
    posterior_from_sufficient_stats_internal(
      data_summ = posterior_stats,
      prior_surv = posterior_prior,
      N_mcmc = 1000,
      single_arm = FALSE
    )
  },
  posterior_from_stats_kernel = {
    set.seed(1001)
    posterior_from_sufficient_stats_kernel_internal(
      data_summ = posterior_stats,
      prior_surv = posterior_prior,
      N_mcmc = 1000,
      single_arm = FALSE
    )
  },
  cox_guarded_auto = {
    cox_wald_test_internal(cox_data)
  },
  cox_public_fallback = {
    cox_wald_test_internal(cox_data, engine = "public")
  },
  impute_success = {
    set.seed(1002)
    impute_data_internal(
      data_in = imputation_data,
      hazard = imputation_hazard,
      end_of_study = end_of_study,
      cutpoints = cutpoints_piecewise,
      type = "success",
      single_arm = FALSE
    )
  },
  impute_futility = {
    set.seed(1003)
    impute_data_internal(
      data_in = imputation_data,
      hazard = imputation_hazard,
      end_of_study = end_of_study,
      cutpoints = cutpoints_piecewise,
      type = "futility",
      single_arm = FALSE
    )
  },
  predictive_imputation_scalar = {
    set.seed(1004)
    scalar_predictive_imputations()
  },
  predictive_imputation_batch = {
    set.seed(1004)
    impute_predictive_draws_internal(
      data_in = imputation_data,
      hazards = imputation_hazards,
      end_of_study = end_of_study,
      cutpoints = cutpoints_piecewise,
      single_arm = FALSE,
      binary_imputation = "event-time",
      check_futility = TRUE
    )
  },
  survival_completed_data_materialized_logrank = {
    materialized_survival_analyses("logrank", 1L)
  },
  survival_completed_data_vectors_logrank = {
    vector_survival_analyses("logrank", 1L)
  },
  survival_completed_data_materialized_cox = {
    materialized_survival_analyses("cox", 1L)
  },
  survival_completed_data_vectors_cox = {
    vector_survival_analyses("cox", 1L)
  },
  survival_completed_data_materialized_bayes = {
    set.seed(1010)
    materialized_survival_analyses("bayes-surv", 300L)
  },
  survival_completed_data_vectors_bayes = {
    set.seed(1010)
    vector_survival_analyses("bayes-surv", 300L)
  },
  binary_completed_data_materialized = {
    materialized_binary_analyses()
  },
  binary_completed_data_counts = {
    count_based_binary_analyses()
  },
  survival_adapt_logrank = {
    set.seed(1005)
    survival_adapt(
      hazard_treatment = trial_hazard_treatment,
      hazard_control = trial_hazard_control,
      cutpoints = trial_cutpoints,
      N_total = 500,
      lambda = 20,
      lambda_time = NULL,
      interim_look = c(250, 375),
      end_of_study = end_of_study,
      prior_surv = c(0.1, 0.1),
      block = 2,
      rand_ratio = c(1, 1),
      prop_loss = 0.2,
      alternative = "less",
      h0 = 0,
      Fn = 0.05,
      Sn = 0.9,
      prob_ha = 0.95,
      N_impute = 100,
      N_mcmc = 100,
      method = "logrank",
      imputed_final = FALSE
    )
  },
  survival_adapt_cox = {
    set.seed(1009)
    survival_adapt(
      hazard_treatment = trial_hazard_treatment,
      hazard_control = trial_hazard_control,
      cutpoints = trial_cutpoints,
      N_total = 500,
      lambda = 20,
      lambda_time = NULL,
      interim_look = c(250, 375),
      end_of_study = end_of_study,
      prior_surv = c(0.1, 0.1),
      block = 2,
      rand_ratio = c(1, 1),
      prop_loss = 0.2,
      alternative = "less",
      h0 = 0,
      Fn = 0.05,
      Sn = 0.9,
      prob_ha = 0.95,
      N_impute = 100,
      N_mcmc = 1,
      method = "cox",
      imputed_final = FALSE
    )
  },
  survival_adapt_bayes = {
    set.seed(1006)
    survival_adapt(
      hazard_treatment = trial_hazard_treatment,
      hazard_control = trial_hazard_control,
      cutpoints = trial_cutpoints,
      N_total = 500,
      lambda = 20,
      lambda_time = NULL,
      interim_look = c(250, 375),
      end_of_study = end_of_study,
      prior_surv = c(0.1, 0.1),
      block = 2,
      rand_ratio = c(1, 1),
      prop_loss = 0.2,
      alternative = "less",
      h0 = 0,
      Fn = 0.05,
      Sn = 0.9,
      prob_ha = 0.95,
      N_impute = 100,
      N_mcmc = 300,
      method = "bayes-surv",
      imputed_final = TRUE
    )
  },
  survival_adapt_bayes_bin = {
    set.seed(1008)
    survival_adapt(
      hazard_treatment = trial_hazard_treatment,
      hazard_control = trial_hazard_control,
      cutpoints = trial_cutpoints,
      N_total = 500,
      lambda = 20,
      lambda_time = NULL,
      interim_look = c(250, 375),
      end_of_study = end_of_study,
      prior_surv = c(0.1, 0.1),
      prior_bin = c(1, 1),
      bin_method = "quadrature",
      block = 2,
      rand_ratio = c(1, 1),
      prop_loss = 0.2,
      alternative = "less",
      h0 = 0,
      Fn = 0.05,
      Sn = 0.9,
      prob_ha = 0.95,
      N_impute = 100,
      N_mcmc = 100,
      method = "bayes-bin",
      imputed_final = TRUE,
      binary_imputation = "bernoulli"
    )
  },
  iterations = iterations,
  check = FALSE,
  filter_gc = FALSE
)

print(benchmark_results, n = Inf)

output_file <- Sys.getenv("GOLDILOCKS_BENCHMARK_OUT", "")
if (nzchar(output_file)) {
  utils::write.csv(
    as.data.frame(benchmark_results),
    output_file,
    row.names = FALSE
  )
  message("Benchmark results written to ", output_file)
}
