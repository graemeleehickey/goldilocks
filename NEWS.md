# goldilocks 0.4.0.9000

## Improvements

* `survival_adapt()` now supports one-sided tests (`alternative = "greater"` or `"less"`) for `method = "cox"` and `method = "logrank"`. The chi-square test remains two-sided only (#20).
* When `method = "chisq"` and `imputed_final = FALSE`, subjects lost to follow-up are now excluded from the final analysis. Previously, LTFU subjects were counted as non-events, which diluted the event rate and biased the chi-square test (#22).
* Removed unused C++ global variables and dead threading code inherited from the deprecated `fastlogranktest` package.
* Replaced all uses of the magrittr pipe (`%>%`) with the base pipe (`|>`) in `summarise_sims()` and `posterior()`, and removed the `dplyr::%>%` re-export from the NAMESPACE.
* `posterior()` now warns when a piecewise interval has zero subjects and data is propagated from an adjacent interval.
* `survival_adapt()` now validates that each `interim_look` in a two-arm design is at least the block size, so an interim look that could enrol a single treatment arm only is rejected as an input error rather than producing an undefined interim posterior.

## Bug fixes

* `posterior()` now propagates data for zero-exposure piecewise intervals *within* each treatment arm. Previously, propagation walked the flat row order of the per-arm summary, so an empty leading interval in the treatment arm (e.g. a look where the treatment arm has no subjects) could copy the control arm's data into the treatment posterior, contaminating the estimate across arms. `posterior()` also now checks up front that the expected treatment arms are present in the supplied data, erroring with an informative message instead of silently returning an all-`NA` posterior slice for an absent arm (which could occur at a small interim look where one arm has no enrolled subjects yet).

* `survival_adapt()` no longer errors when called without interim looks (`interim_look = NULL`). The final analysis previously relied on an undefined loop index variable, which has been replaced with `stage_trial_stopped`.
* `impute_data()` no longer uses hard-coded positional column subsetting (`[, 1:10]`). Temporary columns are now dropped by name, making the function robust to upstream changes in the data frame structure (#26).
* `randomization()` no longer produces `NA` for `next_block` when the loop exhausts all elements of a multi-element `block` vector. The index now wraps around cyclically (#31).
* `analyse_data()` now uses explicit row/column indexing when extracting Cox model results, preventing silent errors if the summary matrix structure changes (#29).
* `analyse_data()` now imports `stats::pnorm`, which is used to compute one-sided `p`-values for the Cox and log-rank tests. Previously this relied on `stats` being attached.
* Removed a dead `loss_to_fu <- NA` assignment in `survival_adapt()` that was shadowed by the `loss_to_fu` column inside `within()` and never used.
* `enrollment()` now correctly selects the enrollment rate at piecewise changepoints. Previously, the rate at exact changepoint boundaries could use the rate from the prior interval (#28).
* `survival_adapt()` no longer adds a systematic perturbation (`sd(time) / 1e4`) to all survival times at interim looks. Instead, only the boundary subject with zero follow-up time is clamped to `.Machine$double.eps` to satisfy `survSplit()` requirements (#24).

## Documentation

* Removed the `est_interim` element from the `survival_adapt()` return-value documentation. This field was documented but never computed or returned.
* Documented the two-stage posterior procedure used when `method = "bayes"` with imputation, clarifying that the imputation model's posterior influences the analysis posterior (#27).
* Clarified `prop_loss` parameter documentation, explaining that LTFU times are drawn from `Uniform(0, t)` and that the event has not yet occurred at the dropout time (#25).

## Housekeeping

* Added unit tests for `enrollment()`, `randomization()`, `pwe_sim()`, `pwe_impute()`, `ppwe()`, `prop_to_haz()`, `sim_comp_data()`, `summarise_sims()`, `analyse_data()`, and `posterior()`.
* Added a package architecture vignette with a function dependency diagram.
* Added `_pkgdown.yml` configuration for a documentation website.
* Added GitHub Actions workflow for pkgdown site deployment.
* Updated GitHub Actions (`actions/checkout`, `actions/upload-artifact`) from v4 to v5 for Node.js 24 compatibility.
* Added `.positai`, `_pkgdown.yml`, and `docs` to `.Rbuildignore` to suppress `R CMD check` NOTEs.
* Removed the unused `appveyor.yml` CI configuration file and its stale `.Rbuildignore` entry.
* Updated README to clarify that the C++ log-rank code was ported from the now-deprecated `fastlogranktest` package.
* Clarified `prior` parameter documentation to explicitly state the Gamma rate parameterization and that the same prior is shared across all piecewise intervals and treatment arms.

# goldilocks 0.4.0

## Main updates

* Because `fastlogranktest` is no longer available on CRAN, a copy of the C++ code and wrapper from the CRAN archive have been included directly into the source code of this package.

## Housekeeping

* Updated GitHub actions workflows
* Updated README badges
* Changed logic check in `summarise_sims.R`

# goldilocks 0.3.1

## Features

* Added a chi-square test for binomial data

## Bugs

* Fixed a small bug in `summarise_sims()`

## Housekeeping

* Restructured code so complete trial data is generated by `sim_comp_data()`
* Restructured code so interim stopping tests are done by `test_stop_success()`
* Restructured code so final analysis is done by `test_final()`
* Added URL links to DESCRIPTION
* Updates to documentation
* Removed appveyor check and badge

# goldilocks 0.3.0

* Added a `NEWS.md` file to track changes to the package.
* First CRAN submission.
