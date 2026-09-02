# Derive piecewise-constant hazard rates from cumulative event probabilities

Converts cumulative event probabilities at specified follow-up times
into the corresponding piecewise-constant hazard rates. This is useful
for expressing data-generation assumptions in terms of clinically
interpretable event probabilities.

## Usage

``` r
prop_to_haz(probs, cutpoints = NULL, endtime)
```

## Arguments

- probs:

  A required numeric vector of finite cumulative event probabilities in
  `[0, 1)` at each cutpoint and at `endtime`, in that order. Its length
  must be one greater than the number of cutpoints. With no cutpoints,
  supply a single probability at `endtime`. Values must be
  non-decreasing and are not recycled.

- cutpoints:

  `NULL` (the default), or a numeric vector of finite, positive,
  strictly increasing interior times at which the event hazard changes.
  `NULL` corresponds to a simple (non-piecewise) exponential model.

- endtime:

  A required single finite, positive numeric value giving the follow-up
  time corresponding to the final element of `probs`. It must be later
  than every cutpoint and use the same time unit.

## Value

A numeric vector of non-negative hazard rates, with one value for each
interval defined by `cutpoints` and `endtime`.

## Details

Given \\J-1\\ interior cutpoints, then there are J intervals defined as:
\\\[s_0, s_1)\\, \\\[s_1, s_2)\\, \\\dots\\, \\\[s\_{J-1}, s\_{J})\\,
with conditions that \\s_0 = 0\\ and \\s_J = \infty\\. Each interval
corresponds to constant hazard \\\lambda_j\\. This is the PWEALL
representation of the continuous generating hazard. Changing the value
at an isolated cutpoint does not alter the cumulative probabilities
calculated here. When observed event times are assigned to analysis
intervals, `goldilocks` uses `(s_{j-1}, s_j]`, matching the survival
counting-process convention, so an event at \\s_j\\ belongs to the
interval ending there.

## Examples

``` r
lambda <- prop_to_haz(0.15, endtime = 36) # 15% probability at 36-months
all.equal(pexp(36, lambda), 0.15)
#> [1] TRUE

# 15% probability at 12-months, and 30% at 24-months
prop_to_haz(c(0.15, 0.30), 12, 24)
#> [1] 0.01354324 0.01617967
PWEALL::pwe(12, prop_to_haz(c(0.15, 0.30), 12, 24), c(0, 12))$dist
#> [1] 0.15
PWEALL::pwe(24, prop_to_haz(c(0.15, 0.30), 12, 24), c(0, 12))$dist
#> [1] 0.3
```
