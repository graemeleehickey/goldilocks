
<!-- README.md is generated from README.Rmd. Please edit that file -->

# goldilocks <img src="man/figures/hex.png" width = "175" height = "200" align="right" />

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/goldilocks)](https://CRAN.R-project.org/package=goldilocks)
[![R build
status](https://github.com/graemeleehickey/goldilocks/workflows/R-CMD-check/badge.svg)](https://github.com/graemeleehickey/goldilocks/actions)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/graemeleehickey/goldilocks?branch=main&svg=true)](https://ci.appveyor.com/project/graemeleehickey/goldilocks)
[![codecov](https://codecov.io/gh/graemeleehickey/goldilocks/branch/main/graph/badge.svg?token=9V6BH1Q4K3)](https://codecov.io/gh/graemeleehickey/goldilocks)
<!-- badges: end -->

The goal of `goldilocks` is to implement the Goldilocks Bayesian
adaptive design proposed by Broglio et al. (2014) for time-to-event
endpoint trials, both one- and two-arm, with an underlying piecewise
exponential hazard model.

The method can be used for a confirmatory trial to select a trial’s
sample size based on accumulating data. During accrual, frequent sample
size selection analyses are made and predictive probabilities are used
to determine whether the current sample size is sufficient or whether
continuing accrual would be futile. The algorithm explicitly accounts
for complete follow-up of all patients before the primary analysis is
conducted. Final analysis tests include the log-rank test, Cox
proportional hazards regression Wald test, and a Bayesian test that
compares the absolute difference in cumulative incidence functions at a
fixed time point.

Broglio et al. (2014) refer to this as a *Goldilocks trial design*, as
it is constantly asking the question, “Is the sample size too big, too
small, or just right?”

## Key benefits

Other software and R packages are available to implement this algorithm.
However, when designing studies it is generally required that many
thousands of trials are simulated to adequately characterize the
operating characteristics, e.g. type I error and power. Hence, a
computationally efficient and fast algorithm is helpful. The
`goldilocks` package takes advantage of many tools to achieve this:

-   Log-rank tests are implemented via the
    [`fastlogranktest`](https://CRAN.R-project.org/package=fastlogranktest)
    package, which uses a lightweight C++ implementation

-   Piecewise exponential simulation is implemented via the
    [`PWEALL`](https://CRAN.R-project.org/package=PWEALL) package, which
    uses a lightweight Fortran implementation

-   Simulation of multiple trials can be performed in parallel using the
    [`pbmcapply`](https://CRAN.R-project.org/package=pbmcapply) package

## References

Broglio KR, Connor JT, Berry SM. Not too big, not too small: a
Goldilocks approach to sample size selection. *Journal of
Biopharmaceutical Statistics*, 2014; **24(3)**: 685–705.

## Installation

You can install the development version of `goldilocks`
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("graemeleehickey/goldilocks")
```
