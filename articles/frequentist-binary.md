# Frequentist binary outcome designs

``` r

library(goldilocks)
```

`goldilocks` provides two frequentist analyses for a two-arm binary
outcome measured at a fixed follow-up time:

- `method = "riskdiff-fm"` uses a Farrington-Manning score test.
- `method = "riskdiff-wald"` retains the plug-in Wald test.

Both methods estimate the treatment-minus-control event-risk difference,

\widehat\Delta = \widehat p\_{\text{treatment}} - \widehat
p\_{\text{control}},

and compare it with the null value or margin supplied through `h0`. All
three alternatives are supported. For example, when events are
undesirable, `alternative = "less"` asks whether the treatment event
risk is lower than the control event risk by enough to cross the chosen
threshold.

## Choosing the analysis

`riskdiff-fm` is the more suitable choice when small samples or rare
outcomes can produce zero events, all events, or little outcome
variation in one arm. Its variance uses maximum likelihood event risks
constrained by the null hypothesis. Under `h0 = 0`, an equal-arm
all-zero or all-one table produces a neutral test result instead of a
zero-variance error. A table with events in only one arm ordinarily
retains a positive constrained variance and is analyzed in the usual
way.

`riskdiff-wald` uses the observed arm risks directly in an unpooled
plug-in variance. It is retained for designs that prespecify that
analysis and for reproducing earlier results. The former
`method = "riskdiff"` name remains a deprecated alias for
`"riskdiff-wald"` and produces a warning. New analyses should prespecify
one of the two procedures explicitly.

The Farrington-Manning calculation has no continuity correction. If its
null-constrained variance is zero, `goldilocks` defines the score
statistic as zero when the observed risk difference equals `h0`,
positive infinity when it is above `h0`, and negative infinity when it
is below `h0`.

## Example design

Suppose the control event probability by 12 months is expected to be
10%, and the treatment event probability is expected to be 5%. The
following small example uses the score test with a one-sided
alternative:

``` r

end_of_study <- 12

result <- survival_adapt(
  hazard_treatment = prop_to_haz(0.05, endtime = end_of_study),
  hazard_control = prop_to_haz(0.10, endtime = end_of_study),
  N_total = 80,
  lambda = 10,
  interim_look = 40,
  end_of_study = end_of_study,
  alternative = "less",
  h0 = 0,
  Fn = 0.05,
  Sn = 0.90,
  prob_ha = 0.975,
  N_impute = 20,
  method = "riskdiff-fm"
)

result
#>   prob_threshold margin alternative N_treatment N_control N_enrolled N_max
#> 1          0.975      0        less          20        20         40    80
#>   post_prob_ha est_final ppp_success stop_futility stop_immediate_success
#> 1    0.9266035      -0.1           0             1                      0
#>   stop_expected_success trial_success stopping_reason decision_time
#> 1                     0         FALSE        futility      4.643332
#>   accrual_stop_time analysis_ready_time planned_completion_time
#> 1          4.643332            16.64333                16.64333
#>   followup_person_time peak_active_followup
#> 1             465.9974                   40
```

For a non-inferiority design, `h0` can be a nonzero risk-difference
margin. If an increase in undesirable events of no more than 5
percentage points is acceptable, use `h0 = 0.05` with
`alternative = "less"`. The score-test variance is then calculated from
event risks constrained to differ by 0.05.

## Pending and missing endpoint outcomes

At an interim look, each pending endpoint is completed from the
piecewise-exponential predictive model before the selected
risk-difference test is applied. `binary_imputation = "event-time"`
samples a conditional event time; `binary_imputation = "bernoulli"`
samples the equivalent fixed-time status directly.

At the final analysis, subjects lost before complete endpoint
ascertainment are excluded when `imputed_final = FALSE`. With
`imputed_final = TRUE`, the package combines risk-difference estimates
and their within-imputation variances using Rubin’s rules. That combined
result is a pooled Wald analysis for either method setting; it is not a
Farrington-Manning test. The same neutral-or-directional convention
prevents an error when its total variance is zero.

The saved results include the selected method and evaluated design
arguments, providing an auditable record of the analysis specification.

## Reference

Farrington CP, Manning G. Test statistics and sample size formulae for
comparative binomial trials with null hypothesis of non-zero risk
difference or non-unity relative risk. *Statistics in Medicine*. 1990;
**9**:1447-1454. <doi:10.1002/sim.4780091208>.
