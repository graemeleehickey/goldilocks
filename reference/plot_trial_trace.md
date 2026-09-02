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

  A required `goldilocks_trial` result from
  [`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md),
  a `goldilocks_interim` result from
  [`evaluate_interim()`](https://graemeleehickey.github.io/goldilocks/reference/evaluate_interim.md),
  or an interim trace data frame.

## Value

The trace data frame, invisibly.
