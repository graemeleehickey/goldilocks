## Submission notes

This is goldilocks 0.6.0, a substantial update following 0.5.0. The main
changes are:

* exact continuous-time enrollment from a piecewise-constant Poisson process,
  with enrollment and hazard schedules now both using internal change points;
* Bayesian beta-binomial and frequentist risk-difference analyses for
  fixed-time binary endpoints;
* separate interim and final piecewise-exponential priors, including
  interval-specific Gamma priors;
* reproducible serial, fork, and PSOCK simulation backends;
* optional interim decision traces and new enrollment, stopping,
  operating-characteristic, and decision-region plots;
* broader input validation, bug fixes, and simulation hot-path performance
  improvements; and
* expanded package, function, and vignette documentation.

The release notes identify the breaking enrollment-schedule and analysis-method
changes, including migration from `lambda_time = 0` to `lambda_time = NULL`,
from `method = "chisq"` to `method = "riskdiff"`, and from `method = "bayes"`
to `method = "bayes-surv"`.

## Current test environment

* local macOS (Tahoe 26.5.2), R 4.5.3, `R CMD check --as-cran`

Before submission, the final release version should also be checked on the
repository's Linux, macOS, and Windows GitHub Actions matrix and on win-builder.

## R CMD check results

0 errors | 0 warnings | 2 notes

* The local check could not verify the current time because the runner could not
  reach an external time server.
* HTML validation and math-rendering checks were skipped locally because a
  sufficiently recent HTML Tidy and the R package `V8` were unavailable.

The package source builds successfully; all 491 testthat expectations pass; all
examples and eight vignettes build and rebuild; the complete pkgdown site
builds; documentation consistency checks pass; and all package URLs validate.

## Reverse dependencies

There are no reverse dependencies on CRAN.
