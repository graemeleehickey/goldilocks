# Single-arm designs with a performance goal

``` r

library(goldilocks)
```

Several other vignettes focus on two-arm randomized designs, although
the Bayesian binary outcome vignette also includes a single-arm example.
Single-arm trials – in which every subject receives the experimental
therapy and the comparator is an external benchmark, often called a
performance goal (PG) or objective performance criterion (OPC) – are
common in early-phase oncology, rare-disease, and proof-of-concept
studies. This vignette shows how to set up a Goldilocks single-arm
design with
[`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md).

Two practical constraints on single-arm designs in this package:

- A single-arm trial is signaled by setting `hazard_control = NULL`.
- `method = "bayes-surv"` supports single-arm survival analyses with
  piecewise-exponential event-time modeling. `method = "bayes-bin"`
  supports single-arm analyses of complete binary outcomes. The
  frequentist tests (`logrank`, `cox`, `chisq`) require two arms and
  will raise an error if used in this mode.

## The decision rule

In a single-arm trial there is no concurrent control, so the “treatment
effect” is replaced by the cumulative event probability on the treatment
arm itself:

$`\text{effect} \;=\; p_{\text{treatment}} \;=\; \Pr(\text{event by end\_of\_study} \mid \text{data}).`$

The argument `h0` plays the role of a benchmark on this scale: a target
failure probability (or, equivalently, $`1 - h_0`$ is a target survival
probability) drawn from external evidence such as a published rate,
registry, or historical cohort. In clinical-trial terminology this
benchmark may be referred to as a performance goal (PG) or objective
performance criterion (OPC). With `alternative = "less"` and `prob_ha`,
the trial declares success when

``` math
\Pr(p_{\text{treatment}} < h_0 \mid \text{data}) \;>\; \texttt{prob\_ha},
```

i.e. when the posterior assigns enough mass to “the experimental therapy
has a lower failure rate than the benchmark”. Choosing
`alternative = "greater"` reverses the direction;
`alternative = "two.sided"` is not allowed for `method = "bayes-surv"`.

The same posterior is used at each interim look to compute the
predictive probability of eventual success, which drives the futility
(`Fn`) and expected-success (`Sn`) stopping rules. Predictive
probabilities are obtained by imputing remaining follow-up from the
posterior predictive distribution of the (piecewise-)exponential model
and re-evaluating the success criterion on each completed dataset.

## Setting up the design

Suppose the existing standard of care has a 30% event probability by 24
months, and we are testing a new agent that we hope will reduce this to
20%. We use an interim look at 50 of 80 enrolled subjects:

``` r

end_of_study <- 24
benchmark <- 0.30                       # external standard-of-care failure rate
target    <- 0.20                       # rate we hope the new therapy achieves

# Convert the target failure rate into a constant hazard (so we can simulate)
ht <- prop_to_haz(probs = target, endtime = end_of_study)
ht
#> [1] 0.009297648
```

Now we run
[`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md):

``` r

out <- survival_adapt(
  hazard_treatment = ht,
  hazard_control   = NULL,              # single-arm
  cutpoints        = 0,
  N_total          = 80,
  lambda           = 5,                 # enrollments per month (constant)
  lambda_time      = 0,
  interim_look     = 50,
  end_of_study     = end_of_study,
  prior            = c(0.1, 0.1),       # Gamma(0.1, 0.1) on the hazard
  prop_loss        = 0.05,
  alternative      = "less",
  h0               = benchmark,         # benchmark failure probability
  Fn               = 0.05,
  Sn               = 0.95,
  prob_ha          = 0.95,
  N_impute         = 50,
  N_mcmc           = 2000,
  method           = "bayes-surv")

out
#>   prob_threshold margin alternative N_treatment N_control N_enrolled N_max
#> 1           0.95    0.3        less          80         0         80    80
#>   post_prob_ha est_final ppp_success stop_futility stop_expected_success
#> 1       0.9985 0.1689997        0.86             0                     0
```

There is no need to supply `block` or `rand_ratio`: they are redundant
in a single-arm design because no randomization is performed.

A few points to highlight in the output:

- `N_control = 0`: no concurrent control was simulated.
- `margin = 0.30`: this is the value of `h0` that the trial is testing
  against. Note that it is on the cumulative-failure scale, not the
  survival scale.
- `est_final` is the posterior mean of $`p_{\text{treatment}}`$ at
  `end_of_study`, *not* a treatment effect relative to control.
- `post_prob_ha` is the posterior probability that
  $`p_{\text{treatment}} < h_0`$.

## Operating characteristics

A single trial does not tell you whether the design is well-calibrated.
To estimate power and type I error, we run the design under each
scenario using
[`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md).
The chunks below are not run when knitting (each takes a few minutes)
but illustrate the workflow:

``` r

# Power: simulate under the alternative (true rate = 0.20)
out_power <- sim_trials(
  N_trials         = 1000,
  hazard_treatment = ht,
  hazard_control   = NULL,
  cutpoints        = 0,
  N_total          = 80,
  lambda           = 5,
  lambda_time      = 0,
  interim_look     = 50,
  end_of_study     = end_of_study,
  prior            = c(0.1, 0.1),
  prop_loss        = 0.05,
  alternative      = "less",
  h0               = benchmark,
  Fn               = 0.05,
  Sn               = 0.95,
  prob_ha          = 0.95,
  N_impute         = 50,
  N_mcmc           = 2000,
  method           = "bayes-surv",
  seed             = 3082)

# Type I error: simulate under the null (true rate = benchmark/PG/OPC = 0.30)
ht_null <- prop_to_haz(probs = benchmark, endtime = end_of_study)
out_t1error <- sim_trials(
  N_trials         = 1000,
  hazard_treatment = ht_null,
  hazard_control   = NULL,
  cutpoints        = 0,
  N_total          = 80,
  lambda           = 5,
  lambda_time      = 0,
  interim_look     = 50,
  end_of_study     = end_of_study,
  prior            = c(0.1, 0.1),
  prop_loss        = 0.05,
  alternative      = "less",
  h0               = benchmark,
  Fn               = 0.05,
  Sn               = 0.95,
  prob_ha          = 0.95,
  N_impute         = 50,
  N_mcmc           = 2000,
  method           = "bayes-surv",
  seed             = 3083)

summarise_sims(list(out_power$sims, out_t1error$sims))
```

Calibration proceeds the same way as for two-arm designs: if the type I
error under the null (where the true rate equals the benchmark) is above
the desired level, raise `prob_ha`; if power is too low, increase
`N_total` or relax the `Fn`/`Sn` thresholds.

## A practical caveat on benchmarks

The validity of a single-arm Goldilocks trial rests entirely on the
benchmark `h0` (the PG or OPC) being a fair representation of the
population the trial is enrolling. Drift in standard of care,
differences in patient mix, and unmeasured confounding all bias the
comparison in a way that randomization would otherwise neutralize. A
Bayesian framework can incorporate uncertainty about the benchmark
itself – e.g. by replacing a fixed `h0` with a prior distribution
informed by historical data – but this is outside the scope of the
simple `h0` scalar that
[`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
exposes, and would require a custom analysis. When in doubt, simulating
the design under several plausible values of the true rate (including
ones near the benchmark) is a useful way to characterize its
sensitivity.

## See also

- The “Two-arm randomized trials” vignette covers the corresponding
  two-arm randomized design with a log-rank decision rule.
- The “Bayesian piecewise-exponential designs” vignette covers the same
  decision rule used here, but in a two-arm setting and with
  non-constant hazards. The piecewise machinery applies directly to
  single-arm trials too (just keep `hazard_control = NULL` and pass a
  per-interval `hazard_treatment` vector).
- [`?survival_adapt`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
  documents all arguments.
