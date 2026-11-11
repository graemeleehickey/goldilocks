## Submission notes

Package: goldilocks 1.1.0. Maintainer: Graeme L. Hickey.

## Correction for the CRAN vignette failure

CRAN reports an error for goldilocks 1.0.0 on r-oldrel-macos-arm64 (R 4.5.2 Patched). Forked R processes terminate with segmentation faults while rebuilding the ADVENT and ANTHEM-HFrEF vignettes; the subsequent missing worker results cause a secondary error when their trial identifiers are collected. Installation, examples, and tests pass on that platform. The check report is available at https://www.R-project.org/nosvn/R.check/r-oldrel-macos-arm64/goldilocks-00check.html.

Evaluated simulations in ADVENT, ANTHEM-HFrEF, and ThermoCool AF now use one core, avoiding forked execution when those vignettes are built or checked. The other evaluated vignette simulations already run sequentially. This change is confined to the vignette execution settings: the package's parallel-processing defaults and error handling are unchanged, and larger examples that are displayed but not evaluated retain their parallel settings.

Simulation sizes, seeds, priors, stopping boundaries, and all other statistical settings are unchanged. In particular, the evaluated scenarios retain 500 trials each in ADVENT, 20 each in ANTHEM-HFrEF, and 50 each in ThermoCool AF. Using one core may increase elapsed build and check time relative to two-core execution, but does not increase the number of simulated trials or predictive calculations. We retain these simulation sizes to preserve the examples' Monte Carlo precision and use sequential execution for portability and reliable vignette rebuilding.

## Verification

GitHub Actions checked commit `62d1670cee30cac1bc156448ef59e7a6fc43f5b2` on 11 November 2026 using released R 4.5.2 on native Apple silicon (`aarch64-apple-darwin20`, macOS Sequoia 15.8.2). The workflow verifies the R version and processor architecture before checking the package. The [successful run and saved check results](https://github.com/graemeleehickey/goldilocks/actions/runs/43057570555) provide the verification record. The original CRAN failure used R 4.5.2 Patched.

The source package was built with `R CMD build --compact-vignettes=gs+qpdf`, including all 18 vignettes. The full `R CMD check --as-cran --timings` completed with **0 errors, 0 warnings, and 0 notes** in 10 minutes 29.1 seconds. All 18 vignettes rebuilt successfully during the check, using 530 seconds of CPU time and 532 seconds elapsed. The retained rebuild log explicitly records successful completion of both `advent.Rmd` and `anthem-hfref.Rmd`. Examples, compiled-code checks, and the PDF and HTML manuals passed. The test suite passed 2,598 expectations with no failures or warnings; eight longer tests were skipped by their existing CRAN conditions. Standard GitHub Actions R setup disables system-clock verification.

The packaged ADVENT and ANTHEM-HFrEF HTML vignettes from the R 4.5.2 build were also inspected visually, including their simulation tables and figures.

Earlier local verification of the one-core correction used R 4.5.3 on macOS Tahoe 26.7.1 (Apple silicon), with standard compiler settings (`R_MAKEVARS_USER=/dev/null`), HTML Tidy 5.8.0, and V8. All 18 vignettes built successfully in 432.83 seconds elapsed. The full check completed with **0 errors, 0 warnings, and 2 notes** in 527.32 seconds elapsed; rebuilding all 18 vignettes used 398 seconds of CPU time and 405 seconds elapsed. The test suite passed 2,500 expectations, with the same eight CRAN skips. The rebuilt ADVENT, ANTHEM-HFrEF, and ThermoCool AF HTML vignettes were also checked visually.

These timings cover the entire 1.1.0 package, including its additional vignettes, and do not isolate the extra elapsed time attributable to switching the three examples to one core.

## Earlier local check notes

* The existing NEJM protocol link in the PROMISE II vignette returned HTTP 403 (Forbidden): https://www.nejm.org/doi/suppl/10.1056/NEJMoa2212754/suppl_file/nejmoa2212754_protocol.pdf. This link is unchanged by the vignette execution correction.
* R reported `unable to verify current time` during the check for future file timestamps.
