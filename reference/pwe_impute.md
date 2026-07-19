# Impute piecewise exponential time-to-event outcomes

Imputation of time-to-event outcomes using the piecewise constant hazard
exponential function conditional on observed exposure.

## Usage

``` r
pwe_impute(time, hazard, cutpoints = NULL, maxtime = NULL)
```

## Arguments

- time:

  vector. The observed time for patient that have had no event or passed
  `maxtime`.

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
