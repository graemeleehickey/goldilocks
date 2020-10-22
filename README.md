
<!-- README.md is generated from README.Rmd. Please edit that file -->

# goldilocks

<!-- badges: start -->

[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/graemeleehickey/goldilocks?branch=master&svg=true)](https://ci.appveyor.com/project/graemeleehickey/goldilocks)
[![Travis build
status](https://travis-ci.com/graemeleehickey/goldilocks.svg?branch=master)](https://travis-ci.com/graemeleehickey/goldilocks)
[![CRAN
status](https://www.r-pkg.org/badges/version/goldilocks)](https://CRAN.R-project.org/package=goldilocks)
[![codecov](https://codecov.io/gh/graemeleehickey/goldilocks/branch/main/graph/badge.svg?token=9V6BH1Q4K3)](https://codecov.io/gh/graemeleehickey/goldilocks)
<!-- badges: end -->

The goal of `goldilocks` is to implement the Golilocks Bayesian adaptive
design proposed by Broglio et al. (2014) for time-to-event endpoint
trials, both one- and two-arm, with an underlying piecewise exponential
hazard model.

The method can be used for a confirmatory trial to select a trial’s
sample size based on accumulating data. During accrual, frequent sample
size selection analyses are made and predictive probabilities are used
to determine whether the current sample size is sufficient or whether
continuing accrual would be futile. The algorithm explicitly accounts
for complete follow-up of all patients before the primary analysis is
conducted.

Broglio et al. (2014) re refer to this as a *Goldilocks trial design*,
as it is constantly asking the question, “Is the sample size too big,
too small, or just right?”

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
