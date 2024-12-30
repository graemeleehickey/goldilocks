## Other notes

Due to the `fastlogranktest` package--a strong dependency for my package--being 
removed from CRAN and the authors not fixing after many months, I have embedded
the C++ code and R wrapper directly into here. I have attributed this to the
original author in the DESCRIPTION file. Note: `fastlogranktest` was distributed
under GPL-3 license.

## Test environments

* local macOS (Sonoma 14.6.1) install, R 4.4.2
* ubuntu (via GitHub actions, release + devel)
* macOS (via GitHub actions, release)
* windows (via appveyor CI, release)
* windows (via GitHub actions, release)
* windows (via win-builder, old + release + devel)

## R CMD check results

0 errors | 0 warnings | 1 notes

Win-Builder NOTE: "checking CRAN incoming feasibility ... NOTE"

## Reverse dependencies

There are no reverse dependencies
