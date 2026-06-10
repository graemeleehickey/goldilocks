## Submission notes

This is a minor release (0.5.0). It contains several bug fixes -- including a
cross-arm contamination bug in the Bayesian posterior at interim looks, a
silent NA-arm posterior when one treatment arm is absent at an interim look,
a missing `stats::pnorm` import, and corrections to `randomization()` and
`enrollment()` -- and adds support for one-sided Cox and log-rank tests. The
local test suite (156 tests) passes, and `R CMD check` is clean apart from
the environmental NOTEs noted below. Two new vignettes have been added.

## Test environments

* local macOS (Tahoe 26.5) install, R 4.5.3
* ubuntu (via GitHub actions, release + devel)
* macOS (via GitHub actions, release)
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
