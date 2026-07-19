# Simulate piecewise exponential time-to-event outcomes

Simulate time-to-event outcomes using the piecewise constant hazard
exponential function.

## Usage

``` r
pwe_sim(n = 1, hazard = 1, cutpoints = NULL, maxtime = NULL)
```

## Arguments

- n:

  integer. The number of random samples to generate. Default is `n = 1`.

- hazard:

  vector. Finite non-negative constant hazard rates for exponential
  failures. If the final rate is zero, `maxtime` must be supplied so
  that subjects without an event can be administratively censored.

- cutpoints:

  finite, positive, strictly increasing vector of interior times at
  which the hazard rate changes. The number of hazard rates must be one
  greater than the number of cutpoints. Use `NULL` for a constant
  hazard.

- maxtime:

  scalar. Optional administrative censoring time. When supplied, it must
  be later than every cutpoint.

## Value

A data frame with simulated follow-up times (`time`) and respective
event indicator (`event`, 1 = event occurred, 0 = censoring).

## Details

See
[`pwe_impute()`](https://graemeleehickey.github.io/goldilocks/reference/pwe_impute.md)
for details.

## Examples

``` r
pwe_sim(10, hazard = c(0.005, 0.001), cutpoints = 3, maxtime = 36)
#>        time event
#> 1  36.00000     0
#> 2  36.00000     0
#> 3  36.00000     0
#> 4  36.00000     0
#> 5  36.00000     0
#> 6  36.00000     0
#> 7  15.53131     1
#> 8  36.00000     0
#> 9  36.00000     0
#> 10 36.00000     0
y <- pwe_sim(n = 1, hazard = c(2.585924e-02, 3.685254e-09),
             cutpoints = 12)
```
