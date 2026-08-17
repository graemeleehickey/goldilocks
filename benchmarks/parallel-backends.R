# Optional maintainer benchmark for sim_trials() parallel backends.
#
# Run from the package root with:
#   source("benchmarks/parallel-backends.R")
#
# PSOCK timings include cluster startup. The benchmark also measures equivalent
# cluster creation and teardown separately, then reports the remaining elapsed
# time as an estimate of compute and dispatch cost. Compare relative results on
# the same macOS or Windows machine.

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

parse_positive_integer_env <- function(name, default) {
  value <- Sys.getenv(name, default)
  parsed <- suppressWarnings(as.integer(strsplit(value, ",", fixed = TRUE)[[
    1
  ]]))
  if (anyNA(parsed) || any(parsed < 1L)) {
    stop(
      name,
      " must contain comma-separated positive integers.",
      call. = FALSE
    )
  }
  unique(parsed)
}

load_goldilocks()

iterations <- parse_positive_integer_env("GOLDILOCKS_PARALLEL_ITERATIONS", "3")
if (length(iterations) != 1L) {
  stop("GOLDILOCKS_PARALLEL_ITERATIONS must be one integer.", call. = FALSE)
}
n_trials_values <- parse_positive_integer_env(
  "GOLDILOCKS_PARALLEL_TRIALS",
  "2,8,32"
)
worker_counts <- sort(unique(c(
  1L,
  parse_positive_integer_env("GOLDILOCKS_PARALLEL_WORKERS", "1,2,4")
)))

end_of_study <- 24
hazard_control <- prop_to_haz(c(0.15, 0.25), 12, end_of_study)
hazard_treatment <- prop_to_haz(c(0.10, 0.18), 12, end_of_study)

workloads <- list(
  small = list(
    hazard_treatment = hazard_treatment,
    hazard_control = hazard_control,
    cutpoints = 12,
    N_total = 60,
    lambda = 10,
    lambda_time = NULL,
    interim_look = NULL,
    end_of_study = end_of_study,
    method = "logrank"
  ),
  logrank = list(
    hazard_treatment = hazard_treatment,
    hazard_control = hazard_control,
    cutpoints = 12,
    N_total = 240,
    lambda = 15,
    lambda_time = NULL,
    interim_look = c(120, 180),
    end_of_study = end_of_study,
    N_impute = 20,
    method = "logrank"
  ),
  bayes_surv = list(
    hazard_treatment = hazard_treatment,
    hazard_control = hazard_control,
    cutpoints = 12,
    N_total = 240,
    lambda = 15,
    lambda_time = NULL,
    interim_look = c(120, 180),
    end_of_study = end_of_study,
    N_impute = 20,
    N_mcmc = 200,
    method = "bayes-surv"
  )
)

backends <- c("auto", "sequential", "psock")
if (.Platform$OS.type != "windows") {
  backends <- c(backends, "fork")
}

cases <- expand.grid(
  workload = names(workloads),
  N_trials = n_trials_values,
  backend = backends,
  ncores = worker_counts,
  stringsAsFactors = FALSE
)
cases <- subset(
  cases,
  (backend == "sequential" & ncores == 1L) |
    (backend != "sequential" & ncores > 1L)
)

benchmark_case <- function(workload, N_trials, backend, ncores) {
  args <- c(
    workloads[[workload]],
    list(
      N_trials = N_trials,
      ncores = ncores,
      backend = backend,
      seed = 4101
    )
  )

  timings <- lapply(seq_len(iterations), function(i) {
    started <- unname(proc.time()[["elapsed"]])
    out <- do.call(sim_trials, args)
    total_seconds <- unname(proc.time()[["elapsed"]]) - started
    metadata <- attr(out, "parallel_metadata", exact = TRUE)
    list(
      total_seconds = total_seconds,
      metadata = metadata
    )
  })
  metadata <- timings[[1]]$metadata
  total_seconds <- vapply(timings, `[[`, numeric(1), "total_seconds")
  setup_seconds <- if (metadata$backend == "psock") {
    vapply(
      seq_len(iterations),
      function(i) {
        started <- unname(proc.time()[["elapsed"]])
        cluster <- parallel::makeCluster(metadata$workers)
        parallel::stopCluster(cluster)
        unname(proc.time()[["elapsed"]]) - started
      },
      numeric(1)
    )
  } else {
    rep.int(0, iterations)
  }
  median_seconds <- stats::median(total_seconds)
  median_setup_seconds <- stats::median(setup_seconds)

  profile_memory <- backend == "sequential"
  memory_bytes <- NA_real_
  if (profile_memory) {
    memory_result <- bench::mark(
      do.call(sim_trials, args),
      iterations = 1L,
      check = FALSE,
      filter_gc = FALSE,
      memory = TRUE
    )
    memory_bytes <- as.numeric(memory_result$mem_alloc)
  }

  data.frame(
    workload = workload,
    N_trials = N_trials,
    backend = backend,
    effective_backend = metadata$backend,
    requested_ncores = ncores,
    effective_workers = metadata$workers,
    selection_reason = metadata$selection_reason,
    median_seconds = median_seconds,
    median_setup_seconds = median_setup_seconds,
    estimated_compute_seconds = max(
      median_seconds - median_setup_seconds,
      0
    ),
    memory_bytes = memory_bytes,
    memory_measurement = if (profile_memory) {
      "parent allocation"
    } else {
      "unavailable for child-process benchmarks"
    },
    stringsAsFactors = FALSE
  )
}

results <- do.call(
  rbind,
  Map(
    benchmark_case,
    cases$workload,
    cases$N_trials,
    cases$backend,
    cases$ncores
  )
)

serial <- results[results$backend == "sequential", ]
names(serial)[names(serial) == "median_seconds"] <- "serial_seconds"
serial <- serial[, c("workload", "N_trials", "serial_seconds")]
results <- merge(results, serial, by = c("workload", "N_trials"), all.x = TRUE)
results$speedup <- results$serial_seconds / results$median_seconds
results$platform <- R.version$platform
results$R_version <- R.version.string
results$os_type <- .Platform$OS.type

results <- results[
  order(
    results$workload,
    results$N_trials,
    results$backend,
    results$requested_ncores
  ),
]
print(results, row.names = FALSE)

output_file <- Sys.getenv("GOLDILOCKS_BENCHMARK_OUT", "")
if (nzchar(output_file)) {
  utils::write.csv(results, output_file, row.names = FALSE)
  message("Benchmark results written to ", output_file)
}
