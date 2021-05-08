## Test environments

* macOS-Catalina (local R installation), R 4.0.3
* ubuntu-20.04 (via GitHub actions, release + devel), R 4.0.3
* windows-latest (via GitHub actions, release)
* macOS-latest (via GitHub actions, release)
* win-builder (release)
* win-builder (devel)
* windows-latest (via appveyor CI, release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Reverse dependencies

There are no reverse dependencies.

## Other notes

This is a resubmission following the CRAN team's comment. Specifically, I have
updated this package per Gregor Seyer's review:

- DESCRIPTION no longer adds extra LICENSE file
- Rd files that are exported all have a return value description.
- The vignette file does not run the code using the 8-cores. The .Rmd file does
evaluate the code due to the computational length of the simulations. Instead,
the simulations have been preprocessed and saved as an .rda file in the
vignettes folder.
