# Simulate exact continuous-time enrollment

Simulate enrollment times from a Poisson process with a
piecewise-constant rate.

## Usage

``` r
enrollment(lambda = 1, N_total, lambda_time = NULL)
```

## Arguments

- lambda:

  finite positive enrollment rates per unit time. Supply one rate for
  each interval defined by `lambda_time`, so `length(lambda)` must equal
  `length(lambda_time) + 1`.

- N_total:

  positive integer total sample size.

- lambda_time:

  `NULL`, or a finite, positive, strictly increasing vector of interior
  times at which the enrollment rate changes. The initial boundary at
  time zero is implicit and must not be supplied. Use `NULL` for a
  constant enrollment rate.

## Value

A non-decreasing numeric vector of `N_total` continuous enrollment
times, measured from first patient in and expressed in the same time
unit used for `lambda_time`. The first value is zero.

## Details

**Major behavior change from goldilocks 0.5.0 and earlier.** Versions
through 0.5.0 generated Poisson counts in unit-time bins. `enrollment()`
returned rebased integer bin times, after which
[`sim_comp_data()`](https://graemeleehickey.github.io/goldilocks/reference/sim_comp_data.md)
added independent uniform jitter and sorted the result. The current
implementation instead generates the continuous arrival times directly
from the exact piecewise-constant Poisson process. Consequently:

- seeded simulations do not reproduce enrollment or downstream trial
  results obtained with version 0.5.0 or earlier;

- enrollment times and operating-characteristic estimates can change,
  particularly when rates are low or a rate change is not an integer
  time;

- the schedule API now contains internal knots only: change
  `lambda_time = 0` to `lambda_time = NULL` for a constant rate, and
  change `lambda_time = c(0, t1, t2)` to `lambda_time = c(t1, t2)` for a
  piecewise rate.

`enrollment()` treats time zero as **first patient in**, matching the
time origin used throughout `goldilocks`. The first returned enrollment
time is therefore fixed at zero. The remaining `N_total - 1` arrivals
form a continuous-time non-homogeneous Poisson process whose rate is
constant between the supplied internal knots. Thus, `lambda` is measured
in enrollments per unit of `lambda_time`; for example, when time is
measured in months, `lambda = 5` means five enrollments per month on
average.

Write the internal knots as \\0 \< \tau_1 \< \cdots \< \tau_K\\, with
\\\tau_0 = 0\\ implicit. The enrollment intensity is

\$\$ \lambda(t) = \lambda_j, \qquad \tau\_{j-1} \le t \< \tau_j, \$\$

for \\j = 1, \ldots, K\\, and \\\lambda(t) = \lambda\_{K+1}\\ after the
final knot. The final rate therefore continues for as long as needed to
reach `N_total`; there is no finite accrual horizon in this function.

Arrivals are generated exactly by the time-rescaling theorem. If \\E_2,
\ldots, E_N\\ are independent unit-rate exponential variables and \\S_i
= \sum\_{k=2}^i E_k\\, then the enrollment times after first patient in
are

\$\$T_1 = 0, \qquad T_i = \Lambda^{-1}(S_i),\$\$

where \\\Lambda(t) = \int_0^t \lambda(u)\\du\\ is the cumulative
enrollment intensity. This construction gives independent Poisson counts
on disjoint calendar intervals, with expected count \\\int_a^b
\lambda(u)\\du\\ over \\(a,b\]\\. For a constant rate, successive
enrollment gaps are independent `Exponential(lambda)` variables.

The rate-change times are measured from first patient in, not from an
earlier site-opening or trial-activation date. If operational delays
before first patient in are important, they must be modelled separately
before using the returned relative times.

For example, `lambda = c(0.3, 0.7, 0.9, 1.2)` with
`lambda_time = c(5, 10, 15)` specifies average enrollment rates of 0.3
over \\\[0,5)\\, 0.7 over \\\[5,10)\\, 0.9 over \\\[10,15)\\, and 1.2
from time 15 onward. Fractional knots such as `lambda_time = 2.5` are
handled exactly; no unit-time binning or post-hoc jitter is used.

## Examples

``` r
# Constant enrollment: first patient at zero, then exponential gaps.
enrollment(lambda = 0.7, N_total = 10)
#>  [1]  0.000000  3.387777  4.343014  4.630902  6.979419 14.372719 15.590170
#>  [8] 16.255547 17.034180 17.068086

# Three internal rate changes define four enrollment intervals.
enrollment(
  lambda = c(0.3, 0.7, 0.9, 1.2),
  N_total = 50,
  lambda_time = c(5, 10, 15)
)
#>  [1]  0.000000  3.248967  6.110215  8.898290 10.595670 12.926809 14.313422
#>  [8] 14.536687 16.280805 16.310522 16.636895 16.951155 19.969558 20.313615
#> [15] 20.417782 20.718054 21.351142 21.800600 22.025401 22.820496 24.306914
#> [22] 24.359537 25.055098 26.824180 27.749159 28.291158 28.948175 29.064916
#> [29] 29.928269 30.088707 31.687406 32.000755 32.673702 33.267533 33.856746
#> [36] 34.857605 34.957459 35.179734 35.347463 35.962991 36.674114 37.483202
#> [43] 38.672304 38.758390 40.229620 43.388860 43.912279 45.002938 46.399839
#> [50] 46.412589

# Fractional change times are supported exactly.
enrollment(
  lambda = c(0.25, 1),
  N_total = 20,
  lambda_time = 2.5
)
#>  [1]  0.000000  1.282874  1.377205  2.890406  3.307968  3.315925  3.564742
#>  [8]  3.832695  4.008059  4.180333  5.485005  6.024266  8.217778  8.544191
#> [15]  8.927070  9.081201  9.456097 10.873318 10.899277 11.086384
```
