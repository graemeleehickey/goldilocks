# Simulate piecewise exponential time-to-event outcomes

Simulates event times from a piecewise-exponential distribution, with
optional administrative censoring at a fixed follow-up time.

## Usage

``` r
pwe_sim(n = 1, hazard = 1, cutpoints = NULL, maxtime = NULL)
```

## Arguments

- n:

  A single non-negative integer giving the number of event times to
  simulate. The default is `1`; `n = 0` returns a zero-row data frame.

- hazard:

  A numeric vector of finite, non-negative event rates, with one value
  per interval defined by `cutpoints`. The default is `1`, giving a
  constant unit rate. If at least one outcome is requested and the final
  rate is zero, `maxtime` must be supplied so that subjects without an
  event can be administratively censored.

- cutpoints:

  `NULL` (the default), or a numeric vector of finite, positive,
  strictly increasing interior times at which the hazard rate changes.
  The number of hazard rates must be one greater than the number of
  cutpoints. Use `NULL` for a constant hazard.

- maxtime:

  `NULL` (the default), or a single finite, positive numeric
  administrative censoring time. When supplied, it must be later than
  every cutpoint.

## Value

A data frame with one row per simulated subject and columns `time`, the
event or censoring time, and `event`, coded `1` for an event and `0` for
administrative censoring.

## Details

PWEALL represents the generating hazard with pieces closed on the left
and open on the right. Because the event-time distribution is
continuous, the value of the hazard at an isolated cutpoint does not
alter the cumulative hazard, distribution, or generated samples. When
realized event times are later assigned to analysis intervals,
`goldilocks` follows the survival counting-process convention, open on
the left and closed on the right, so an event exactly at a cutpoint
belongs to the interval ending there. See
[`pwe_impute()`](https://graemeleehickey.github.io/goldilocks/reference/pwe_impute.md)
for the conditional sampling details.

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
