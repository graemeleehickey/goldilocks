## Submission notes

Package: goldilocks 1.1.0. Maintainer: Graeme L. Hickey.

## Correction for the CRAN vignette failure

CRAN reports an error for goldilocks 1.0.0 on r-oldrel-macos-arm64 (R 4.5.2 Patched). Forked R processes terminate with segmentation faults while rebuilding the ADVENT and ANTHEM-HFrEF vignettes; the subsequent missing worker results cause a secondary error when their trial identifiers are collected. Installation, examples, and tests pass on that platform. The check report is available at https://www.R-project.org/nosvn/R.check/r-oldrel-macos-arm64/goldilocks-00check.html.

Evaluated simulations in ADVENT, ANTHEM-HFrEF, and ThermoCool AF now use one core, avoiding forked execution when those vignettes are built or checked. The other evaluated vignette simulations already run sequentially. This change is confined to the vignette execution settings: the package's parallel-processing defaults and error handling are unchanged, and larger examples that are displayed but not evaluated retain their parallel settings.

Simulation sizes, seeds, priors, stopping boundaries, and all other statistical settings are unchanged. In particular, the evaluated scenarios retain 500 trials each in ADVENT, 20 each in ANTHEM-HFrEF, and 50 each in ThermoCool AF. Using one core may increase elapsed build and check time relative to two-core execution, but does not increase the number of simulated trials or predictive calculations. We retain these simulation sizes to preserve the examples' Monte Carlo precision and use sequential execution for portability and reliable vignette rebuilding.

## Verification

Verification used R 4.5.3 on macOS Tahoe 26.7.1 (Apple silicon), with standard compiler settings (`R_MAKEVARS_USER=/dev/null`), HTML Tidy 5.8.0, and V8. R 4.5.2, used by the failing CRAN check, and R-devel are not installed locally.

The source package was built with `R CMD build --compact-vignettes=gs+qpdf`; all 18 vignettes built successfully in 432.83 seconds elapsed. The rebuilt ADVENT, ANTHEM-HFrEF, and ThermoCool AF HTML vignettes were also checked visually.

The full `R CMD check --as-cran --timings` completed with **0 errors, 0 warnings, and 2 notes** in 527.32 seconds elapsed. All 18 vignettes rebuilt successfully; this step used 398 seconds of CPU time and 405 seconds elapsed. Examples, compiled-code checks, and the PDF and HTML manuals passed. The test suite passed 2,500 expectations with no failures or warnings; eight longer tests were skipped by their existing CRAN conditions.

These timings cover the entire 1.1.0 package, including its additional vignettes, and do not isolate the extra elapsed time attributable to switching the three examples to one core.

## Check notes

* The existing NEJM protocol link in the PROMISE II vignette returned HTTP 403 (Forbidden): https://www.nejm.org/doi/suppl/10.1056/NEJMoa2212754/suppl_file/nejmoa2212754_protocol.pdf. This link is unchanged by the vignette execution correction.
* R reported `unable to verify current time` during the check for future file timestamps.
