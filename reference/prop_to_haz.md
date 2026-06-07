# Estimate plausible piecewise constant hazard rates from summary summary event proportions

Given estimates of the event probability at one or more fixed times, the
corresponding piecewise hazard rates can be determined through
closed-form formulae. This utility function can be useful when
simulating trial datasets with plausible event rates.

## Usage

``` r
prop_to_haz(probs, cutpoints = 0, endtime)
```

## Arguments

- probs:

  vector. Probabilities of the event (i.e. cumulative incidence
  probabilities) at one or more time point. If only a single value is
  given, then it is assumed that this is the probability at the
  `endtime`.

- cutpoints:

  vector. Times at which the baseline hazard changes. Default is
  `cutpoints = 0`, which corresponds to a simple (non-piecewise)
  exponential model.

- endtime:

  scalar. Time at which final element in `probs` corresponds to.
  Typically this would be the study endpoint time.

## Value

Vector of constant hazard rates for each time piece defined by
`cutpoints`.

## Details

Given \\J-1\\ internal cut-points, then there are J intervals defined
as: \\\[s_0, s_1)\\, \\\[s_1, s_2)\\, \\\dots\\, \\\[s\_{J-1},
s\_{J})\\, with conditions that \\s_0 = 0\\ and \\s_J = \infty\\. Each
interval corresponds to constant hazard \\\lambda_j\\.

## Examples

``` r
lambda <- prop_to_haz(0.15, endtime = 36) # 15% probability at 36-months
all.equal(pexp(36, lambda), 0.15)
#> [1] TRUE

# 15% probability at 12-months, and 30% at 24-months
prop_to_haz(c(0.15, 0.30), c(0, 12), 24)
#> [1] 0.01354324 0.01617967
PWEALL::pwe(12, prop_to_haz(c(0.15, 0.30), c(0, 12), 24), c(0, 12))$dist
#> [1] 0.15
PWEALL::pwe(24, prop_to_haz(c(0.15, 0.30), c(0, 12), 24), c(0, 12))$dist
#> [1] 0.3
```
