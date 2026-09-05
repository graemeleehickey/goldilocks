# Statistical workflow and reproducibility

## Overview

This vignette describes how the statistical components of a `goldilocks`
design fit together. Its purpose is to support protocol development,
simulation reports, and independent review: the assumptions used to
generate trial data, the models used for interim prediction, the
decision thresholds, and the final analysis should be distinguishable
and prespecified.

The package supports three related analyses:

1.  [`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md)
    simulates repeated trials to estimate operating characteristics;
2.  [`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
    simulates and evaluates one adaptive trial; and
3.  [`evaluate_interim()`](https://graemeleehickey.github.io/goldilocks/reference/evaluate_interim.md)
    applies the prespecified interim rule to an observed trial data cut.

The same posterior-predictive calculation and ordered stopping rule are
used for simulated and observed interim analyses.

## From assumptions to operating characteristics

### Data-generating assumptions

For simulation,
[`sim_comp_data()`](https://graemeleehickey.github.io/goldilocks/reference/sim_comp_data.md)
generates enrollment times, randomized treatment assignments, event
times, and loss to follow-up. The first participant enrolls at time
zero; subsequent enrollment follows a piecewise-constant Poisson
process. Event times follow arm-specific piecewise-exponential
distributions.

The event-time distribution used to simulate trials need not equal the
model used for interim prediction. `generation_cutpoints` defines the
data-generating hazard intervals, whereas `cutpoints` defines the
intervals used for posterior estimation, predictive imputation, and
Bayesian survival analysis. This separation permits sensitivity analyses
for model misspecification without changing the prespecified analysis
model.

### Interim prediction

At look \ell, let n\_\ell be the number enrolled and let N\_{\max} be
the maximum sample size. The observed interim data update the Gamma
priors for the piecewise-exponential hazards. Posterior-predictive
simulation then estimates

P\_{n\_\ell} = \Pr(\text{final success after follow-up of the current
cohort} \mid \mathcal{D}^{\mathrm{obs}}\_\ell)

and

P\_{\max,\ell} = \Pr(\text{final success after enrollment to } N\_{\max}
\mid \mathcal{D}^{\mathrm{obs}}\_\ell).

Each predictive replicate completes pending outcomes under a posterior
draw of the event-time hazards and applies the prespecified
completed-data analysis. The proportions of successful replicates
estimate P\_{n\_\ell} and P\_{\max,\ell}.

The ordered decision rule is

d\_\ell = \begin{cases} \text{declare immediate success}, &
\widehat{P}\_{n\_\ell} \> Q\_\ell, \\ \text{stop accrual for expected
success}, & S\_\ell \< \widehat{P}\_{n\_\ell} \le Q\_\ell, \\
\text{declare binding futility}, & \widehat{P}\_{n\_\ell} \le S\_\ell
\text{ and } \widehat{P}\_{\max,\ell} \< F\_\ell, \\ \text{continue
enrollment}, & \text{otherwise.} \end{cases}

The package requires Q\_\ell \ge S\_\ell. With the default Q\_\ell = 1,
immediate success is disabled. All comparisons are strict, so equality
with a boundary does not cross it.

`N_impute` controls the number of posterior-predictive replicates. For
Bayesian completed-data analyses, `N_mcmc` controls the posterior draws
within each replicate. The reported Monte Carlo standard errors and
exact binomial bounds describe numerical uncertainty in the predictive
probabilities; the interim decision itself uses the point estimate.

For fixed-horizon binary analyses, event counts and denominators by
treatment arm are sufficient statistics. Carrying these sufficient
statistics into the completed-data analysis gives the same
risk-difference or beta-binomial analysis as participant-level endpoint
records.

## Bayesian survival calculation

Let d\_{aj} and y\_{aj} denote the observed event count and person-time
in arm a and interval j. With the independent prior

\lambda\_{aj} \sim \operatorname{Gamma}(\alpha\_{0aj},\beta\_{0aj}),

Gamma-exponential conjugacy gives

\lambda\_{aj} \mid \mathcal{D}^{\mathrm{obs}}\_\ell \sim
\operatorname{Gamma}( \alpha\_{0aj}+d\_{aj}, \beta\_{0aj}+y\_{aj} ).

For analysis cut-points 0\<c_1\<\cdots\<c\_{J-1}\<\tau, define interval
widths

\mathbf{w} = (c_1,c_2-c_1,\ldots,\tau-c\_{J-1}).

A posterior hazard draw implies cumulative hazard and event probability

H_a(\tau)=\sum\_{j=1}^{J}\lambda\_{aj}w_j, \qquad
p_a(\tau)=1-\exp\\-H_a(\tau)\\.

These widths span the fixed endpoint horizon `end_of_study`; they are
not shortened to the longest follow-up observed at an interim look. If
no participant has contributed information to a later interval,
`empty_interval = "prior"` leaves its posterior equal to its prior. This
makes prior-predictive assessment particularly important when little
late follow-up is expected at early looks.

## Completed-data analysis

The completed-data method determines whether each predictive replicate,
and when required the final trial data, meets the success criterion.

| `method` | Estimand and analysis | With imputation at the final analysis | Without final imputation |
|:---|:---|:---|:---|
| `bayes-surv` | Posterior treatment-minus-control event probability, or the treatment event probability in a single-arm design | Average posterior summaries across completed imputations | Analyze observed right-censored follow-up |
| `bayes-bin` | Beta-binomial posterior for fixed-horizon event status | Average posterior summaries across completed imputations | Exclude participants without complete endpoint ascertainment |
| `cox` | Log hazard ratio from a Cox model | Pool estimates and variances using Rubin’s rules | Analyze observed right-censored follow-up |
| `riskdiff-wald` | Treatment-minus-control event-risk difference using a Wald test | Pool estimates and variances using Rubin’s rules | Exclude participants without complete endpoint ascertainment |
| `riskdiff-fm` | Treatment-minus-control event-risk difference using a Farrington-Manning score test | Pool estimates and variances using Rubin’s rules, yielding a pooled Wald analysis | Exclude participants without complete endpoint ascertainment |
| `logrank` | Difference between survival distributions using a log-rank test | Not available because an imputation-pooling rule has not been specified | Analyze observed right-censored follow-up |

For a frequentist method, the success measure is 1-p; for a Bayesian
method, it is the posterior probability of the prespecified alternative.
In both cases, success requires the measure to be strictly greater than
`prob_ha`.

The imputation model and completed-data analysis model are deliberately
distinct. For example, `method = "bayes-bin"` uses the
piecewise-exponential model to impute endpoint status for participants
whose endpoint is pending, then applies a beta-binomial model to the
completed binary outcomes. Both sets of assumptions should therefore be
examined in sensitivity analyses.

## Summaries for design evaluation

For one simulated trial or an observed data cut, the interim decision
history shows the two predictive probabilities, their Monte Carlo
uncertainty, the three decision thresholds, and the action at each
completed look. For repeated simulations,
[`summarise_sims()`](https://graemeleehickey.github.io/goldilocks/reference/summarise_sims.md)
estimates power or type I error, probabilities of each stopping outcome,
sample-size summaries, and Monte Carlo uncertainty.
[`summarise_calendar_time()`](https://graemeleehickey.github.io/goldilocks/reference/summarise_calendar_time.md)
adds trial duration, accrual duration, analysis readiness, and follow-up
burden.

The planning functions
[`prop_to_haz()`](https://graemeleehickey.github.io/goldilocks/reference/prop_to_haz.md)
and
[`ppwe()`](https://graemeleehickey.github.io/goldilocks/reference/ppwe.md)
connect clinically interpretable event probabilities to
piecewise-exponential hazards.
[`plot_enrollment()`](https://graemeleehickey.github.io/goldilocks/reference/plot_enrollment.md)
displays the expected enrollment trajectory and planned interim
milestones. These quantities are useful checks that the numerical design
corresponds to the assumptions intended for the protocol.

## Prespecification and reproducibility

Several conventions should be stated explicitly in a simulation report
or statistical analysis plan:

- `treatment = 0` denotes control and `treatment = 1` denotes treatment;
- named arm-specific inputs use `control`, then `treatment`;
- `generation_cutpoints` governs event-time generation, whereas
  `cutpoints` governs prediction and analysis;
- analysis intervals use the survival counting-process convention
  (\text{start},\text{stop}\] when assigning observed events;
- `Qn`, `Sn`, `Fn`, `prob_ha`, `N_impute`, and `N_mcmc` are part of the
  prespecified decision algorithm; and
- a recorded simulation seed permits exact reproduction of the Monte
  Carlo study across the supported computing options.

Operating characteristics should be evaluated over clinically plausible
null and alternative scenarios, including nuisance parameters that may
affect the amount of information available at interim looks. Numerical
Monte Carlo error should be reported alongside every estimated operating
characteristic.
