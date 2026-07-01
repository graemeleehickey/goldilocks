# Optional maintainer benchmarks for simulation hot paths.
#
# Run from the package root with:
#   source("benchmarks/hot-paths.R")

if (!requireNamespace("bench", quietly = TRUE)) {
  stop("Install the 'bench' package to run these benchmarks.", call. = FALSE)
}

load_goldilocks <- function() {
  if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    library(goldilocks)
  }
}

load_goldilocks()

ppwe_internal <- getFromNamespace("ppwe", "goldilocks")
haz_to_prop_internal <- getFromNamespace("haz_to_prop", "goldilocks")
posterior_internal <- getFromNamespace("posterior", "goldilocks")
impute_data_internal <- getFromNamespace("impute_data", "goldilocks")

iterations <- as.integer(Sys.getenv("GOLDILOCKS_BENCHMARK_ITERATIONS", "5"))
if (is.na(iterations) || iterations < 1L) {
  stop("GOLDILOCKS_BENCHMARK_ITERATIONS must be a positive integer.", call. = FALSE)
}

set.seed(4242)

cutpoints_piecewise <- c(0, 6, 12)
end_of_study <- 36

hazard_matrix <- matrix(
  stats::rgamma(15000, shape = 2, rate = 80),
  ncol = length(cutpoints_piecewise)
)

posterior_draws <- array(
  stats::rgamma(30000, shape = 2, rate = 80),
  dim = c(5000, length(cutpoints_piecewise), 2)
)

posterior_data <- data.frame(
  time = pmin(stats::rexp(600, rate = 0.03), end_of_study),
  event = stats::rbinom(600, size = 1, prob = 0.65),
  treatment = rep(0:1, each = 300)
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
  dim = c(1, length(cutpoints_piecewise), 2)
)

trial_cutpoints <- c(0, 4, 8)
trial_hazard_control <- prop_to_haz(c(0.20, 0.30, 0.35), trial_cutpoints, end_of_study)
trial_hazard_treatment <- prop_to_haz(c(0.10, 0.18, 0.22), trial_cutpoints, end_of_study)

benchmark_results <- bench::mark(
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
  posterior_piecewise = {
    set.seed(1001)
    posterior_internal(
      data = posterior_data,
      cutpoints = cutpoints_piecewise,
      prior = c(0.1, 0.1),
      N_mcmc = 1000,
      single_arm = FALSE
    )
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
  survival_adapt_logrank = {
    set.seed(1004)
    survival_adapt(
      hazard_treatment = trial_hazard_treatment,
      hazard_control = trial_hazard_control,
      cutpoints = trial_cutpoints,
      N_total = 500,
      lambda = 20,
      lambda_time = 0,
      interim_look = c(250, 375),
      end_of_study = end_of_study,
      prior = c(0.1, 0.1),
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
      imputed_final = TRUE
    )
  },
  survival_adapt_bayes = {
    set.seed(1005)
    survival_adapt(
      hazard_treatment = trial_hazard_treatment,
      hazard_control = trial_hazard_control,
      cutpoints = trial_cutpoints,
      N_total = 500,
      lambda = 20,
      lambda_time = 0,
      interim_look = c(250, 375),
      end_of_study = end_of_study,
      prior = c(0.1, 0.1),
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
      method = "bayes",
      imputed_final = TRUE
    )
  },
  iterations = iterations,
  check = FALSE,
  filter_gc = FALSE
)

print(benchmark_results)

output_file <- Sys.getenv("GOLDILOCKS_BENCHMARK_OUT", "")
if (nzchar(output_file)) {
  utils::write.csv(as.data.frame(benchmark_results), output_file, row.names = FALSE)
  message("Benchmark results written to ", output_file)
}
