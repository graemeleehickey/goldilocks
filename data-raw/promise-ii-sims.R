# Reproduce the PROMISE II vignette results from the repository root.
# Author: Graeme L. Hickey
# Usage: Rscript data-raw/promise-ii-sims.R [trials=10000] [draws=10000] [cores=4]
# Intermediate scenario files permit resuming the same calculation after interruption.

pkgload::load_all(".", quiet = TRUE)
source("inst/examples/promise-ii.R")
args <- commandArgs(trailingOnly = TRUE)
N_trials <- if (length(args) >= 1L) as.integer(args[1]) else 10000L
N_impute <- if (length(args) >= 2L) as.integer(args[2]) else 10000L
ncores <- if (length(args) >= 3L) as.integer(args[3]) else 4L
helper_hash <- unname(tools::md5sum("inst/examples/promise-ii.R"))
cache_dir <- file.path("cache", "promise-ii")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

scenarios <- data.frame(
  name = c(paste0("primary_", c(54, 65:70, 75)), "slower_null", "faster_null",
           "fixed_null", "fixed_target"),
  afs = c(.54, .65, .66, .67, .68, .69, .70, .75, .54, .54, .54, .70),
  accrual = c(rep("expected", 8), "slower", "faster", "expected", "expected"),
  enrollment_model = c(rep("poisson", 10), "fixed", "fixed"),
  seed = 11012026L + seq_len(12)
)
simulations <- vector("list", nrow(scenarios))
names(simulations) <- scenarios$name
for (i in seq_len(nrow(scenarios))) {
  scenario <- scenarios[i, ]
  path <- file.path(cache_dir, paste0(scenario$name, ".rds"))
  request <- list(scenario = scenario, N_trials = N_trials, N_impute = N_impute,
                  helper_hash = helper_hash)
  cached <- if (file.exists(path)) readRDS(path) else NULL
  if (!is.null(cached) && identical(cached$request, request)) {
    simulations[[i]] <- cached$result
    message("Reusing ", scenario$name)
  } else {
    message("Running ", scenario$name, ": ", N_trials, " trials, ", N_impute, " draws")
    simulations[[i]] <- promise_ii_simulate(
      afs = scenario$afs, N_trials = N_trials, N_impute = N_impute,
      seed = scenario$seed, accrual = scenario$accrual,
      enrollment_model = scenario$enrollment_model, ncores = ncores
    )
    saveRDS(list(request = request, result = simulations[[i]]), path)
  }
  print(promise_ii_summary(simulations[[i]]))
}
results <- list(
  simulations = simulations, scenarios = scenarios,
  reference = promise_ii_reference(),
  metadata = list(author = "Graeme L. Hickey", helper_md5 = helper_hash,
                  generator_md5 = unname(tools::md5sum("data-raw/promise-ii-sims.R")),
                  R_version = R.version.string,
                  package_version = as.character(packageVersion("goldilocks")),
                  source = "SAP Rev. 5, Tables 3 and 5; planned maximum N = 120",
                  created = format(Sys.time(), tz = "UTC", usetz = TRUE),
                  session = capture.output(sessionInfo()))
)
saveRDS(results, "vignettes/promise-ii-results.rds", compress = "xz")
# The retained results include every trial; the per-scenario copies are disposable.
unlink(file.path(cache_dir, paste0(scenarios$name, ".rds")))
if (!length(list.files(cache_dir, all.files = TRUE, no.. = TRUE))) unlink(cache_dir, recursive = TRUE)
