# Package architecture

## Overview

The `goldilocks` package has three ways into the adaptive calculation:

1.  [`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md)
    simulates many independent trials and delegates each one to
    [`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md);
2.  [`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
    generates and evaluates one complete simulated trial; and
3.  [`evaluate_interim()`](https://graemeleehickey.github.io/goldilocks/reference/evaluate_interim.md)
    applies the same interim decision calculation to an externally
    observed data cut.

All three routes ultimately use the same interim decision engine. This
vignette separates the package-level flow from the more detailed
Bayesian-survival calculation so that the performance paths do not
obscure the statistical model.

In the diagrams, **blue** nodes are exported functions, **grey** nodes
are internal orchestration or checked analysis functions, and **gold**
nodes are trusted internal kernels. A trusted kernel receives canonical
objects that have already been validated and normalized by `goldilocks`;
it does not define a different statistical analysis.

## Trial simulation and interim evaluation

[`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md)
first resolves sequential, fork, or PSOCK execution and, when a seed is
supplied, constructs an independent `"L'Ecuyer-CMRG"` stream for every
requested trial. Trial-level errors are isolated and retained in the
returned failure table. Each successful task calls
[`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
with the same normalized design but its own random-number stream.

[`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
calls
[`sim_comp_data()`](https://graemeleehickey.github.io/goldilocks/reference/sim_comp_data.md)
to generate enrollment, treatment assignments, event times, and loss to
follow-up. `generation_cutpoints` defines only the event-time
data-generating partition. The separate `cutpoints` argument defines
posterior estimation, predictive imputation, and Bayesian-survival
analysis.

At every simulated look,
[`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
prepares the observed data cut and calls `evaluate_interim_core()`.
[`evaluate_interim()`](https://graemeleehickey.github.io/goldilocks/reference/evaluate_interim.md)
instead validates an external data cut, derives potential future
allocation from `N_total` and `rand_ratio`, and then calls the same core
function. Consequently, simulated and external interim decisions share
posterior prediction, imputation, completed-data analysis, threshold
comparisons, diagnostics, and trace construction.

The core draws `N_impute` hazards from the observed-data posterior and
passes the complete array to `impute_predictive_draws()`. That internal
batch calculation generates subject-by-draw time and event matrices for
currently pending follow-up and, when futility is enabled, future
subjects through the maximum sample size. It stores only rows requiring
imputation. `test_stop_success()` then materializes and analyzes one
matrix column at a time, so the completed-data analysis models and
decision rules remain unchanged. The scalar `impute_data()` helper
remains in active use for final-analysis imputation and as the reference
implementation.

Within an interim batch, random inputs are generated in posterior-draw
order; for each draw the order is current treatment, current control,
future treatment, then future control. All predictive imputations are
generated before completed-data analysis begins. This remains exactly
reproducible under [`set.seed()`](https://rdrr.io/r/base/Random.html)
and across supported serial and parallel backends. Because Bayesian
completed-data analyses previously consumed random numbers between
imputations, the same seed can produce different Bayesian results than
in earlier package versions. The resulting success indicators estimate
the predictive probabilities used for expected-success and futility
decisions. Exact Monte Carlo bounds are retained as diagnostics;
decisions use the point estimates.

## Bayesian-survival checked and trusted paths

Bayesian-survival analysis uses the same Gamma piecewise-exponential
model on both paths below. The difference is where validation occurs and
which intermediate objects are materialized.

The checked route accepts the package’s supported prior specifications
and defensively verifies sufficient-statistic columns, arm-interval
combinations, counts, and empty intervals while normalizing arbitrary
valid row order. It remains the route used by ordinary `analyse_data()`
calls, including observed-data Bayesian-survival final analyses.
`haz_to_prop()` returns treatment, control, and effect draws and uses
the exported
[`ppwe()`](https://graemeleehickey.github.io/goldilocks/reference/ppwe.md)
calculation for piecewise hazards.

The trusted route is used only after `goldilocks` has produced canonical
arm-by-interval sufficient statistics and normalized the Gamma priors.
It retains inexpensive invariant checks, draws from exactly the same
Gamma posterior, and consumes random numbers in exactly the same order
as the checked route. It then calculates only the treatment-effect draws
needed by the decision rule, avoiding repeated cutpoint validation and
temporary probability data frames.

For analysis cutpoints 0 \< c_1 \< \cdots \< c\_{J-1} \< \tau, where
\tau is `end_of_study`, `endpoint_interval_widths()` returns

\mathbf{w} = (c_1,\\ c_2-c_1,\\ \ldots,\\ \tau-c\_{J-1}).

For arm a, a posterior hazard draw is converted to cumulative hazard and
event probability using

H_a(\tau) = \sum\_{j=1}^{J}\lambda\_{aj}w_j, \qquad p_a(\tau) = 1 -
\exp\\-H_a(\tau)\\.

These widths describe the fixed analysis horizon, not the maximum
follow-up currently observed. Observed events and person-time enter
separately through the sufficient statistics used in the Gamma
posterior. If nobody has entered a late interval,
`empty_interval = "prior"` leaves that interval prior-driven;
`"propagate"` applies the configured adjacent-interval rule; and
`"error"` stops. Shortening \mathbf{w} to observed follow-up would
change the estimand rather than correct the cumulative-hazard
calculation.

## Completed-data analysis methods

`analyse_data()` dispatches a completed dataset according to `method`.
Interim prediction tests each completed dataset independently before
averaging its binary success indicators.

| `method` | Completed-data analysis | `imputed_final = TRUE` | `imputed_final = FALSE` |
|:---|:---|:---|:---|
| `bayes-surv` | Gamma piecewise-exponential posterior; effect is the treatment-minus-control event probability, or the treatment probability in a single-arm design | Analyze each completed imputation through the trusted sufficient-statistics and effect kernels, then average posterior summaries | Analyze observed right-censored data through the checked path |
| `bayes-bin` | Beta-binomial posterior for completed binary endpoint status | Analyze each completed imputation and average posterior summaries | Exclude subjects lost before complete endpoint ascertainment |
| `cox` | Cox Wald test for the log hazard ratio | Pool estimates and within-imputation variances using Rubin’s rules | Analyze observed right-censored data |
| `riskdiff` | Wald test for the treatment-minus-control event-risk difference | Pool estimates and within-imputation variances using Rubin’s rules | Exclude subjects lost before complete endpoint ascertainment |
| `logrank` | Log-rank test | Not supported because no imputation-pooling rule is implemented | Analyze observed right-censored data |

The posterior-imputation model and completed-data analysis model are
deliberately modular. In particular, `method = "bayes-bin"` still uses
the piecewise-exponential Gamma model to impute pending outcomes, then
applies the separate beta-binomial model to completed endpoint statuses.

## Results and post-processing

## Exported utilities and retained outputs

The exported low-level utilities are:

- [`enrollment()`](https://graemeleehickey.github.io/goldilocks/reference/enrollment.md)
  for exact continuous-time accrual;
- [`randomization()`](https://graemeleehickey.github.io/goldilocks/reference/randomization.md)
  for blocked two-arm allocation;
- [`pwe_sim()`](https://graemeleehickey.github.io/goldilocks/reference/pwe_sim.md)
  and
  [`pwe_impute()`](https://graemeleehickey.github.io/goldilocks/reference/pwe_impute.md)
  for unconditional and conditional piecewise-exponential event times;
- [`ppwe()`](https://graemeleehickey.github.io/goldilocks/reference/ppwe.md)
  for piecewise-exponential cumulative event probabilities;
- [`prop_to_haz()`](https://graemeleehickey.github.io/goldilocks/reference/prop_to_haz.md)
  for converting cumulative event probabilities to piecewise hazards;
  and
- [`sim_comp_data()`](https://graemeleehickey.github.io/goldilocks/reference/sim_comp_data.md)
  for generating a complete simulated trial dataset.

The three main entry points retain enough information for audit and
downstream displays:

- [`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
  returns the one-row trial summary and optionally a compact interim
  trace, with normalized arguments, decision design, prior diagnostics,
  enrollment design, and calendar metrics attached;
- [`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md)
  returns successful trial summaries, isolated failure records, optional
  traces, RNG and parallel metadata, and the common retained design; and
- [`evaluate_interim()`](https://graemeleehickey.github.io/goldilocks/reference/evaluate_interim.md)
  returns the decision, predictive probabilities, Monte Carlo
  diagnostics, allocation and data-cut diagnostics, posterior
  diagnostics, and a compatible one-look trace.

Post-processing functions consume these retained results without
rerunning the trial:

- [`summarise_sims()`](https://graemeleehickey.github.io/goldilocks/reference/summarise_sims.md)
  calculates operating characteristics and Monte Carlo uncertainty;
- [`summarise_calendar_time()`](https://graemeleehickey.github.io/goldilocks/reference/summarise_calendar_time.md)
  summarizes trial duration, accrual, analysis readiness, follow-up, and
  interim timing;
- [`summarise_trial_trace()`](https://graemeleehickey.github.io/goldilocks/reference/summarise_trial_trace.md)
  condenses a single trial’s stopping path;
- [`plot_enrollment()`](https://graemeleehickey.github.io/goldilocks/reference/plot_enrollment.md)
  displays expected or retained accrual and analysis milestones;
- [`plot_trial_trace()`](https://graemeleehickey.github.io/goldilocks/reference/plot_trial_trace.md)
  displays one trial’s predictive-probability path;
- [`plot_sim_stopping()`](https://graemeleehickey.github.io/goldilocks/reference/plot_sim_stopping.md)
  summarizes stopping outcomes and sample sizes;
- [`plot_sim_decisions()`](https://graemeleehickey.github.io/goldilocks/reference/plot_sim_decisions.md)
  displays interim decision regions across simulated trials; and
- [`plot_sim_ocs()`](https://graemeleehickey.github.io/goldilocks/reference/plot_sim_ocs.md)
  compares operating characteristics across scenarios.

## Architectural invariants

Several conventions link the layers above:

- data use `treatment = 0` for control and `treatment = 1` for
  treatment;
- normalized arm-specific vectors and priors use canonical `control`,
  `treatment` order, while posterior hazard arrays retain treatment in
  slice 1 and control in slice 2 for compatibility;
- `generation_cutpoints` governs event-time generation only, whereas
  `cutpoints` governs posterior estimation, predictive imputation,
  analysis priors, and endpoint probability calculations;
- analysis intervals follow the survival counting-process convention
  `(start, stop]` when observed events are assigned to sufficient
  statistics;
- outer predictive decisions use `N_impute`, completed-data Bayesian
  Monte Carlo uses `N_mcmc`, and exact bounds are diagnostic rather than
  decision-changing; and
- interim RNG order is posterior hazards, the complete
  predictive-imputation batch in draw/current/future/arm order, and then
  completed-data analyses; and
- trusted kernels may bypass repeated rich-input validation only after
  package-owned callers have constructed canonical inputs. General entry
  points retain checked behavior.
