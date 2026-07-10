# Plot predictive probabilities and enrollment at interim looks

Draws three base-R panels showing the predictive probability of success
if accrual stops now, the predictive probability of success at the
maximum sample size, and enrollment and observed events by treatment
arm. Thresholds and early stopping decisions are marked on the
probability panels.

## Usage

``` r
plot_trial_trace(x)
```

## Arguments

- x:

  A goldilocks_trial object or an interim trace data frame.

## Value

The trace data frame, invisibly.
