## Other notes

Due to the `fastlogranktest` package--a strong dependency for my package--being 
removed from CRAN and the authors not fixing after many months, I have embedded
the C++ code and R wrapper directly into here. I have attributed this to the
original author in the DESCRIPTION file. Note: `fastlogranktest` was distributed
under GPL-3 license.

## Test environments

* local macOS (Tahoe 26.5) install, R 4.5.3
* ubuntu (via GitHub actions, release + devel)
* macOS (via GitHub actions, release)
* windows (via appveyor CI, release)
* windows (via GitHub actions, release)
* windows (via win-builder, old + release + devel)

## R CMD check results

0 errors | 0 warnings | 1 notes

Win-Builder NOTE: "checking CRAN incoming feasibility ... NOTE"

Local NOTE: "checking for future file timestamps ... NOTE: unable to verify
current time". This is an environment artifact caused by the check being unable
to reach an external time server during the local run; it is unrelated to the
package contents.

## Reverse dependencies

There are no reverse dependencies
