# Run from the package root: Rscript benchmarks/rmst.R
# Optional argument: number of benchmark iterations (default 10).
pkgload::load_all(".", quiet = TRUE)
iterations <- as.integer(commandArgs(trailingOnly = TRUE)[1])
if (is.na(iterations)) {
  iterations <- 10L
}
set.seed(1420)
completed <- sim_comp_data(
  hazard_control = 0.1,
  hazard_treatment = 0.06,
  N_total = 200,
  end_of_study = 12
)
censored <- completed
censored$time[1:10] <- pmin(censored$time[1:10], 2)
censored$event[1:10] <- 0
estimate <- function(data, engine) {
  rmst_estimate(data$time, data$event, data$treatment, 8, engine)
}
stopifnot(isTRUE(all.equal(
  estimate(completed, "auto"),
  estimate(completed, "survfit"),
  tolerance = 1e-12
)))
print(bench::mark(
  completed_auto = estimate(completed, "auto"),
  completed_km = estimate(completed, "survfit"),
  censored_km = estimate(censored, "survfit"),
  iterations = iterations,
  check = FALSE
)[, c("expression", "median", "mem_alloc")])
trial <- function(method) {
  set.seed(1421)
  survival_adapt(
    hazard_control = 0.1,
    hazard_treatment = 0.06,
    N_total = 200,
    end_of_study = 12,
    rmst_tau = 8,
    interim_look = c(100, 150),
    N_impute = 100,
    method = method,
    alternative = if (method == "rmst") "greater" else "less"
  )
}
print(bench::mark(
  rmst_trial = trial("rmst"),
  cox_trial = trial("cox"),
  iterations = iterations,
  check = FALSE
)[, c("expression", "median", "mem_alloc")])
