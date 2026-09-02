# Calculate endpoint event probabilities from piecewise hazards

Calculates the cumulative event probability at a fixed follow-up time
for one or more sets of piecewise-constant hazard rates.

## Usage

``` r
ppwe(hazard, end_of_study, cutpoints = NULL)
```

## Arguments

- hazard:

  A required numeric matrix of finite, non-negative hazard rates. Rows
  represent parameter sets, such as posterior draws, and columns
  represent the intervals defined by `cutpoints`. The number of columns
  must equal `length(cutpoints) + 1`, and at least one row is required.

- end_of_study:

  A required single finite, positive numeric time at which the
  cumulative event probability is evaluated. It must be greater than
  every cutpoint.

- cutpoints:

  `NULL` (the default), or a numeric vector of finite, positive,
  strictly increasing interior times at which the hazard rate changes.
  The number of hazard rates must be one greater than the number of
  cutpoints. Use `NULL` for a constant hazard.

## Value

A numeric vector of event probabilities in `[0, 1]`, with one value for
each row of `hazard`.

## Details

The cumulative probability depends on interval durations, so the value
assigned to an isolated cutpoint has no effect. PWEALL represents its
generating hazard with pieces closed on the left and open on the right.
When `goldilocks` assigns realized event times to analysis intervals, it
instead uses the survival counting-process convention, open on the left
and closed on the right.
