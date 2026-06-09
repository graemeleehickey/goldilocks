# Simulate piecewise exponential time-to-event outcomes

Simulate time-to-event outcomes using the piecewise constant hazard
exponential function.

## Usage

``` r
pwe_sim(n = 1, hazard = 1, cutpoints = 0, maxtime = NULL)
```

## Arguments

- n:

  integer. The number of random samples to generate. Default is `n=1`.

- hazard:

  vector. The constant hazard rates for exponential failures.

- cutpoints:

  vector. The change-point vector indicating time when the hazard rates
  change. Note the first element of `cutpoints` should always be 0.

- maxtime:

  scalar. Maximum time before end of study.

## Value

A data frame with simulated follow-up times (`time`) and respective
event indicator (`event`, 1 = event occurred, 0 = censoring).

## Details

See
[`pwe_impute`](https://graemeleehickey.github.io/goldilocks/reference/pwe_impute.md)
for details.

## Examples

``` r
pwe_sim(10, hazard = c(0.005, 0.001), cutpoints = c(0, 3), maxtime = 36)
#>        time event
#> 1  36.00000     0
#> 2  36.00000     0
#> 3  19.86483     1
#> 4  36.00000     0
#> 5  36.00000     0
#> 6  36.00000     0
#> 7  36.00000     0
#> 8  36.00000     0
#> 9  36.00000     0
#> 10 36.00000     0
y <- pwe_sim(n = 1, hazard = c(2.585924e-02, 3.685254e-09),
             cutpoints = c(0, 12))
```
