# Bayesian piecewise-exponential designs

``` r

library(goldilocks)
```

The “Two-arm randomized trials” vignette uses `method = "logrank"` and a
single constant hazard per arm. That is convenient when the
proportional-hazards assumption is reasonable and the underlying event
rate is well-approximated by an exponential distribution. In practice,
neither assumption always holds. This vignette shows how to set up a
Goldilocks design with:

- a **piecewise-exponential** hazard, in which the constant hazard rate
  changes at one or more cut-points; and
- a **Bayesian decision rule** (`method = "bayes"`) using a posterior
  probability threshold on the cumulative-failure-probability scale.

We use a small example so that the simulation can be run while reading
the vignette; the live
[`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
chunk takes around 10–20 seconds and is cached on first knit.

## When piecewise hazards help

Survival in many settings – post-surgical mortality, transplant
outcomes, oncology with early treatment toxicity – shows a higher hazard
early and a lower steady-state hazard later. A single exponential rate
cannot capture both regimes; using the average rate over-states early
survival and under-states late survival. A piecewise-exponential model
with a single internal cut-point is often a good compromise, fitting one
hazard for the early interval and another for everything afterwards.

The package handles piecewise hazards through two related arguments:

- `cutpoints`: the times at which the hazard changes. The first element
  must be `0` and represents the start of the first interval; subsequent
  elements are the *internal* cut-points. With `J` intervals there are
  `J` elements in `cutpoints` and `J` hazard pieces per arm. (So
  `cutpoints = c(0, 6)` produces *two* intervals, `[0, 6)` and
  `[6, \infty)`, in contrast to the default `cutpoints = 0` used in the
  “Two-arm randomized trials” vignette, which produces a single interval
  and a non-piecewise exponential.)
- `hazard_treatment` and `hazard_control`: vectors of length `J` giving
  the constant hazard for each interval.

The helper
[`prop_to_haz()`](https://graemeleehickey.github.io/goldilocks/reference/prop_to_haz.md)
maps cumulative event probabilities at given times into the
corresponding piecewise hazards.

## Setting up the design

Suppose the primary endpoint is overall survival at 24 months. Based on
prior data, we expect the **control** arm to have:

- 30% mortality by 6 months (an “early high-risk” window), and
- 50% mortality by 24 months overall.

The treatment is hypothesised to attenuate the early-period hazard,
leaving the late-period hazard unchanged. Concretely, we target:

- 18% mortality by 6 months in the treatment arm, and
- 40% mortality by 24 months overall.

We set the cut-point at 6 months and translate these proportions into
piecewise hazards:

``` r

cutpoints <- c(0, 6)         # one internal cut at 6 months -> two intervals
end_of_study <- 24

hc <- prop_to_haz(probs = c(0.30, 0.50), cutpoints = cutpoints, endtime = end_of_study)
ht <- prop_to_haz(probs = c(0.18, 0.40), cutpoints = cutpoints, endtime = end_of_study)

round(rbind(control = hc, treatment = ht), 4)
#>             [,1]   [,2]
#> control   0.0594 0.0187
#> treatment 0.0331 0.0174
```

The first column is the hazard during `[0, 6)` months and the second
column is the hazard from 6 months onward. Both arms have a higher
hazard in the early window. We can sanity-check that the implied
survival probabilities at 24 months match what we specified by running
the cumulative incidence computation back through
[`ppwe()`](https://graemeleehickey.github.io/goldilocks/reference/ppwe.md):

``` r

ppwe(hazard       = matrix(hc, nrow = 1),
     cutpoints    = cutpoints,
     end_of_study = end_of_study)
#> [1] 0.5
ppwe(hazard       = matrix(ht, nrow = 1),
     cutpoints    = cutpoints,
     end_of_study = end_of_study)
#> [1] 0.4
```

These should be 0.50 and 0.40 respectively (modulo rounding).

## The Bayesian decision rule

With `method = "bayes"`,
[`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
puts independent Gamma$`(\alpha, \beta)`$ priors on each piecewise
hazard rate (one per interval, per arm) and combines them with the
observed exposure time and event counts to obtain a closed-form Gamma
posterior on each $`\lambda_j`$. Posterior draws of $`\lambda_j`$ are
pushed through the piecewise-exponential cumulative incidence function
to obtain posterior draws of the cumulative-failure probability $`p`$ at
`end_of_study` for each arm.

The decision rule is one-sided. The treatment effect is defined as

``` math
\Delta = p_{\text{treatment}} - p_{\text{control}},
```

i.e. the difference in failure (not survival) probabilities at
`end_of_study`. On the survival scale this is equivalent to
$`S_{\text{control}} - S_{\text{treatment}}`$, with a negative
$`\Delta`$ corresponding to the treatment having higher survival. With
`alternative = "less"` and a margin `h0` (default `0`), the trial
declares success at the final analysis when

$`\Pr(\Delta < h_0 \mid \text{data}) \;>\; \texttt{prob\_ha}.`$

Because a beneficial treatment has a *lower* failure probability,
`alternative = "less"` is the appropriate choice here.
(`method = "bayes"` does not allow `alternative = "two.sided"` – it
raises an error.)

The same posterior is also used at each interim look to compute the
predictive probability of eventual success. Imputed completions are
drawn from the posterior predictive distribution of the
piecewise-exponential model for subjects still under follow-up, and the
analysis is repeated on each imputed dataset. The fraction of
imputations that would declare success at the maximum sample size is the
predictive probability of success, which drives the futility and
expected-success stopping rules through `Fn` and `Sn`.

We use a weakly informative Gamma$`(0.1, 0.1)`$ prior on every hazard
component:

``` r

prior <- c(0.1, 0.1)   # shape and rate of the Gamma prior on each lambda_j
```

## A single simulated trial

We will run one trial under the alternative hypothesis to illustrate the
mechanics. We choose `interim_look = 60` (the minimum required is
`max(block) = 4`) and a constant accrual rate of 5 enrolments per month
(`lambda_time = 0` keeps the rate constant; piecewise accrual is also
supported).

``` r

out <- survival_adapt(
  hazard_treatment = ht,
  hazard_control   = hc,
  cutpoints        = cutpoints,
  N_total          = 100,
  lambda           = 5,                # enrolments per month
  lambda_time      = 0,                # constant accrual rate
  interim_look     = 60,
  end_of_study     = end_of_study,
  prior            = prior,
  block            = 4,
  rand_ratio       = c(1, 1),
  prop_loss        = 0.05,
  alternative      = "less",
  h0               = 0,
  Fn               = 0.05,
  Sn               = 0.95,
  prob_ha          = 0.975,
  N_impute         = 50,
  N_mcmc           = 2000,
  method           = "bayes")

out
#>   prob_threshold margin alternative N_treatment N_control N_enrolled N_max
#> 1          0.975      0        less          50        50        100   100
#>   post_prob_ha   est_final ppp_success stop_futility stop_expected_success
#> 1        0.556 -0.01536969         0.1             0                     0
```

The output reports the posterior probability of the alternative at the
final (or stopped) analysis (`post_prob_ha`), the posterior mean
treatment effect on the cumulative-failure scale (`est_final`), the
predictive probability of success (`ppp_success`), and indicators for
whether the trial stopped early for futility or expected success.

## Notes on the piecewise model

Two practical considerations are worth flagging:

1.  **Empty intervals at interim looks.** Early interim looks may have
    no subjects with follow-up reaching the later piecewise intervals.
    The package handles this by propagating exposure time and event
    counts from the nearest non-empty interval *within the same arm*,
    and emits a warning when it does so. This supplies a placeholder
    Gamma posterior for an interval with no information; it is a
    fallback when an interim look has not yet generated follow-up in
    later intervals, not a substantive estimate of those intervals’
    hazards. By the final analysis, all intervals will typically be
    populated.

2.  **Number of cut-points.** Each additional cut-point adds two hazard
    parameters to estimate (one per arm). With limited interim data this
    can make individual interval posteriors diffuse. In our experience,
    one or two well-motivated cut-points (e.g. tied to a clinical
    milestone) is usually sufficient; finer partitions tend to add
    variance without commensurate bias reduction.

## Sensitivity to the cut-point specification

If you suspect a piecewise structure but are unsure where the cut-point
should sit, a useful sensitivity check is to fix the data-generating
hazards (the truth) and vary the *analysis* cut-points – that is, what
[`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
is told to model. The simulator generates trial data from
`hazard_treatment` and `hazard_control` evaluated against the supplied
`cutpoints`, so to compare analysis-model choices on a common ground you
would run separate
[`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md)
simulations under the same data-generating truth, varying the analytic
cut-points each time. If the operating characteristics are similar
across choices, the design is robust to the cut-point specification. If
they differ markedly, the cut-point becomes a design decision worth
justifying in the protocol.

A simpler – but very different – comparison is the equivalent design
that is *both* simulated and analysed under a single constant hazard per
arm matched to the overall 24-month proportions:

``` r

hc_flat <- prop_to_haz(0.50, endtime = end_of_study)   # control, single hazard
ht_flat <- prop_to_haz(0.40, endtime = end_of_study)   # treatment, single hazard

out_flat <- survival_adapt(
  hazard_treatment = ht_flat,
  hazard_control   = hc_flat,
  cutpoints        = 0,
  N_total          = 100,
  lambda           = 5,
  lambda_time      = 0,
  interim_look     = 60,
  end_of_study     = end_of_study,
  prior            = prior,
  block            = 4,
  rand_ratio       = c(1, 1),
  prop_loss        = 0.05,
  alternative      = "less",
  h0               = 0,
  Fn               = 0.05,
  Sn               = 0.95,
  prob_ha          = 0.975,
  N_impute         = 50,
  N_mcmc           = 2000,
  method           = "bayes")
```

Note that this changes both the simulated data-generating process *and*
the analysis model, so any difference in operating characteristics
conflates the two effects. It is most useful when the question is “how
would the trial behave if the world really were a single exponential?”
rather than “how robust is my analysis cut-point?”.

## See also

- The “Two-arm randomized trials” vignette covers the same design
  machinery using a single exponential hazard and a log-rank decision
  rule.
- [`?survival_adapt`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
  documents all arguments, including the requirement that each
  `interim_look` in a two-arm design be at least the block size.
- [`?prop_to_haz`](https://graemeleehickey.github.io/goldilocks/reference/prop_to_haz.md)
  and
  [`?ppwe`](https://graemeleehickey.github.io/goldilocks/reference/ppwe.md)
  document the conversion between event proportions and piecewise
  hazards.
