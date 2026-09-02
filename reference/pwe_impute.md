# Impute piecewise exponential time-to-event outcomes

Draws an event time from a piecewise-exponential distribution
conditional on a subject remaining event-free through the observed
follow-up time, with optional administrative censoring.

## Usage

``` r
pwe_impute(time, hazard, cutpoints = NULL, maxtime = NULL)
```

## Arguments

- time:

  A required numeric vector of finite, non-negative event-free follow-up
  times for subjects who have not had an event. When `maxtime` is
  supplied, no value may exceed it. A zero-length vector returns a
  zero-row data frame. Values are not recycled against other arguments.

- hazard:

  A required numeric vector of finite, non-negative event rates, with
  one value per interval defined by `cutpoints`. If the final rate is
  zero, `maxtime` must be supplied.

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

A data frame with one row per subject and columns `time`, the imputed
event or censoring time, and `event`, coded `1` for an event and `0` for
administrative censoring.

## Details

If a subject is event-free at time \\s \< t\\, then the conditional
probability is

\$\$F\_{T \| s}(t \| s) = P(T \le t \| T \> s) = \frac{F(t) - F(s)}{1 -
F(s)}\$\$

where \\F(\cdot)\\ is the cumulative distribution function of the
piecewise exponential (PWE) distribution. Equivalently, \\F(t) = 1 -
S(t)\\, where `S(t)` is the survival function. If \\U \sim Unif(0, 1)\\,
then we can generate an event time (conditional on being event free up
until \\s\\) as

\$\$F^{-1}(U(1 - F(s)) + F(s))\$\$

If \\s = 0\\, this is equivalent to a direct unconditional sample from
the PWE distribution.

PWEALL represents the generating hazard with pieces closed on the left
and open on the right. Its cumulative distribution is continuous at
every cutpoint, so this endpoint choice does not affect imputation. For
assigning realized event times to analysis intervals, `goldilocks` uses
the survival counting-process convention, open on the left and closed on
the right; an event exactly at a cutpoint belongs to the interval ending
there.

## Examples

``` r
pwe_impute(time = c(3, 4, 5), hazard = c(0.002, 0.01), cutpoints = 12)
#>        time event
#> 1  37.24678     1
#> 2 247.16279     1
#> 3 181.45567     1
pwe_impute(time = c(3, 4, 5), hazard = c(0.002, 0.01), cutpoints = 12,
           maxtime = 36)
#>       time event
#> 1 36.00000     0
#> 2 35.24618     1
#> 3 36.00000     0
pwe_impute(time = 19.621870008, hazard = c(2.585924e-02, 3.685254e-09),
           cutpoints = 12, maxtime = 36)
#>   time event
#> 1   36     0
```
