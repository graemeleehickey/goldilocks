# Technical details of prespecified endpoint horizons

The endpoint-horizon functions allow the endpoint time to change
according to a prespecified schedule, conditional on the trial reaching
each interim look. The success calculation uses the endpoint assigned to
that look. The futility calculation predicts success at the eventual
maximum-sample-size endpoint. When these endpoints differ, their
predictive models use separate follow-up data cuts.

This vignette describes the statistical calculations and their
implementation in
[`sim_horizon_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_horizon_trials.md),
[`sim_horizon_paths()`](https://graemeleehickey.github.io/goldilocks/reference/sim_horizon_paths.md),
and
[`apply_horizon_rules()`](https://graemeleehickey.github.io/goldilocks/reference/apply_horizon_rules.md).
[Extending an endpoint from 12 to 18
months](https://graemeleehickey.github.io/goldilocks/articles/endpoint-horizons-example.md)
provides a worked example. The [technical details of the Goldilocks
design](https://graemeleehickey.github.io/goldilocks/articles/technical-methods.md)
cover the underlying fixed-horizon methods.

## Design specification

Let the interim enrollment sizes be n_1 \< \cdots \< n_L \< N\_{\max},
with success endpoint times H_1 \le \cdots \le H_L \le H_F. Here, H_F is
the endpoint for the maximum-N final analysis and for every futility
prediction. All participants have planned outcome collection through H_C
\ge H_F from enrollment, including participants recruited before the
endpoint changes.

| Quantity | R argument | Interpretation |
|:---|:---|:---|
| n\_\ell, H\_\ell | `endpoint_schedule` | Data frame with columns `N` and `horizon`; one row per interim |
| N\_{\max} | `N_total` | Maximum enrollment; not an additional interim look |
| H_F | `final_horizon` | Maximum-N endpoint and futility prediction horizon |
| H_C | `followup_horizon` | Planned collection duration for each participant |
| Supported endpoint times | `endpoint_horizons` | Horizons evaluated and retained by [`sim_horizon_paths()`](https://graemeleehickey.github.io/goldilocks/reference/sim_horizon_paths.md) |
| Dropout reference time | `dropout_horizon` | Time at which `prop_loss` specifies the dropout probability |

The schedule is applied before evaluating each reached look. It does not
use an interim estimate or predictive probability to choose whether to
extend the endpoint. A terminal stopping decision prevents later
scheduled changes from taking effect. If the trial continues to maximum
enrollment, its final analysis uses H_F after endpoint ascertainment;
reaching N\_{\max} does not itself trigger an interim analysis.

## Two data cuts at the same calendar date

Let E_i be enrollment time, T_i the event time from enrollment, and D_i
the independent dropout time. The calendar date of look \ell is t\_\ell
= E\_{n\_\ell}. For an enrolled participant and an analysis horizon H,
the available data are

X\_{i\ell}(H) = \min\\T_i, D_i, H, t\_\ell-E_i\\, \qquad
\delta\_{i\ell}(H) = I\\T_i \le \min(D_i,H,t\_\ell-E_i)\\.

The success data set uses H=H\_\ell; the futility data set uses H=H_F.
Participants not yet enrolled contribute no observed events or exposure
to either model. Their outcomes are generated only for the maximum-N
prediction. The implementation gives a participant enrolled exactly at
the look a negligible positive exposure to accommodate the survival
calculations.

For example, an event at month 15 is represented as event-free
administrative censoring at month 12 in a 12-month success calculation.
It can enter an 18-month futility calculation only if month 15 has
already been reached at that participant’s current follow-up. A dropout
at month 15 does not make the 12-month endpoint missing. A dropout at
month 8 leaves the outcome after month 8 unobserved in both
calculations.

**Figure 1. Separate success and futility predictions.** When the
current and final horizons are equal, the data cuts coincide and the
predictions share one fitted hazard model. When they differ, fitting one
model to the longer data cut and changing only the endpoint of the
completed-data analysis would allow late follow-up to inform early
success; the implementation therefore fits each data cut separately.

## Predictive hazard models and completed-data analyses

### Gamma updating within each data cut

The analysis model has piecewise-constant hazards \lambda\_{zj} for arm
z and interval j, with independent shape-rate priors

\lambda\_{zj} \sim \operatorname{Gamma}(a\_{zj}, b\_{zj}).

For a given look and horizon H, let d\_{zj\ell}(H) and y\_{zj\ell}(H) be
the event count and total exposure in the interval after applying both
the calendar and endpoint cuts. With `empty_interval = "prior"`, the
predictive posterior is

\lambda\_{zj}\mid\mathcal D\_\ell(H) \sim
\operatorname{Gamma}\\a\_{zj}+d\_{zj\ell}(H),\\
b\_{zj}+y\_{zj\ell}(H)\\.

`cutpoints` defines the master analysis partition. At horizon H, only
cutpoints strictly less than H are retained, and the last interval ends
at H. The corresponding leading columns of the interval-specific priors
are retained in the same order. A cutpoint equal to H does not create an
interval beyond the endpoint. Later intervals and their priors are not
combined with earlier ones.

An interval with no observed exposure remains governed by its own prior
under the default policy. This can matter substantially for the
longer-horizon futility calculation at early looks. The other available
policies are `"propagate"`, the existing adjacent-interval heuristic,
and `"error"`. [Bayesian piecewise-exponential
designs](https://graemeleehickey.github.io/goldilocks/articles/bayes-piecewise.md)
discusses predictive and analysis priors in more detail.

### Predictive completion and success scores

For each of the M= `N_impute` posterior draws, the success calculation
imputes outstanding event times in the enrolled cohort through H\_\ell.
The futility calculation imputes outstanding event times among enrolled
participants and generates outcomes for the remaining participants up to
N\_{\max}, all through H_F. Already observed events remain fixed.
Imputation uses event times, including for binary analyses;
independently imputed Bernoulli outcomes at different horizons are not
used.

Each completed cohort is analyzed using the prespecified `method`. A
Bayesian analysis returns the posterior probability of its alternative;
a frequentist analysis returns 1-p for the specified test. Write these
scores as s^{\mathrm{now}}\_{\ell m} and s^{\max}\_{\ell m}. For c=
`prob_ha`, the predictive probabilities are

\widehat P\_\ell = \frac{1}{M}\sum\_{m=1}^{M} I(s^{\mathrm{now}}\_{\ell
m}\>c), \qquad \widehat P\_{\max,\ell} = \frac{1}{M}\sum\_{m=1}^{M}
I(s^{\max}\_{\ell m}\>c).

`prior_surv` supplies the predictive hazard priors. For
`method = "bayes-surv"`, each completed-data analysis forms a fresh
posterior using `prior_surv_final` and the completed cohort. It does not
use the predictive posterior as a new prior, which would count the
observed information twice. For `method = "bayes-bin"`, the completed
endpoint indicators are analyzed with `prior_bin`; the Gamma model
supplies the imputation. `prior_surv_final` also governs hazard-based
imputation at the actual final analysis when requested.

The optional `horizon_settings` data frame specifies horizon-specific
`h0`, `rmst_tau`, or `prior_bin`. Without overrides, the common null
value and beta prior are used at every horizon, while the RMST
restriction equals the endpoint horizon. This allows, for example,
different single-arm performance goals at different endpoint times. The
analysis method and alternative direction remain common across horizons.

The success and futility probabilities integrate over their respective
posteriors. Draw number m from one data cut is not a shared latent event
trajectory with draw number m from another data cut. In contrast, the
underlying simulated trial histories used to evaluate eventual outcomes
are shared across horizons.

## Stopping and final analysis

At look \ell, the thresholds are Q\_\ell= `Qn`, S\_\ell= `Sn`, and
F\_\ell= `Fn`, with Q\_\ell\ge S\_\ell. The first applicable rule is
used.

| Condition | Decision | Endpoint and subsequent analysis |
|:---|:---|:---|
| \widehat P\_\ell\>Q\_\ell | Immediate success | Declare success at H\_\ell without a required final analysis |
| S\_\ell\<\widehat P\_\ell\le Q\_\ell | Expected success | Stop accrual; analyze the enrolled cohort at H\_\ell after endpoint ascertainment |
| \widehat P\_\ell\le S\_\ell and \widehat P\_{\max,\ell}\<F\_\ell | Binding futility | Stop and classify the adaptive trial as unsuccessful |
| Otherwise | Continue | Reach the next scheduled look, or enroll to maximum N and analyze at H_F |

All comparisons are strict. `Qn = 1` disables immediate success,
`Sn = 1` disables expected-success stopping, and `Fn = 0` or `NULL`
disables futility. The Monte Carlo bounds in a trace are diagnostics and
do not replace the point estimates in these comparisons.

**Figure 2. The stopping sequence for a prespecified schedule.**
Futility is evaluated only when enabled and neither success rule
applies. The diagram describes the decisions reconstructed by
[`apply_horizon_rules()`](https://graemeleehickey.github.io/goldilocks/reference/apply_horizon_rules.md)
as well as direct simulation. Complete-path simulation evaluates all
planned looks regardless of these hypothetical decisions.

An expected-success stop at H\_\ell fixes the endpoint for that cohort.
Failure at its subsequent final analysis does not reopen the trial at a
later endpoint. The futility probability concerns the final analysis at
N\_{\max} and H_F; it is not the probability of declaring success at any
remaining interim under the entire adaptive policy.

Final analysis with dropout follows the existing method-specific choices
for `imputed_final`. With `FALSE`, survival analyses retain
right-censoring, while binary analyses exclude participants without
complete endpoint status. With `TRUE`, supported methods impute missing
endpoint outcomes. Independent dropout does not by itself make
complete-case binary analysis unbiased. These choices and the
restrictions on log-rank and Farrington-Manning final imputation are
described in the [technical methods
vignette](https://graemeleehickey.github.io/goldilocks/articles/technical-methods.html#loss-to-follow-up-at-the-final-analysis).

## Shared histories, dropout, and timing

Event times and dropout times are generated once per simulated trial
through `followup_horizon`. Final cohort data are re-censored at each
endpoint; changing the endpoint does not redraw either process.
`generation_cutpoints` and the arm-specific generating hazards must
describe the entire collection period, and need not match the analysis
partition.

If the arm-specific dropout probability is q_z by reference time H_D,
then

D_i\mid Z_i=z \sim \operatorname{Exponential}(\eta_z),\qquad
\eta_z=-\log(1-q_z)/H_D.

Here, H_D= `dropout_horizon`. This hazard remains fixed across endpoint
times. The probability of dropout by the collection horizon is
1-\exp(-\eta_z H_C); the proportion actually censored by dropout can be
smaller because events may occur first.

Trial results distinguish the endpoint attached to the stopping cohort
from continued collection through H_C:

| Selected-endpoint measure | Continued-collection measure | Meaning |
|:---|:---|:---|
| `analysis_ready_time` | `collection_ready_time` | Latest observed event or censoring availability, no earlier than accrual stop |
| `planned_completion_time` | `planned_collection_time` | Accrual-stop time plus the corresponding planned horizon |
| `followup_person_time` | `collection_person_time` | Total observed event-free follow-up through the corresponding horizon |

The collection measures concern the enrolled cohort, including after an
early stop. They do not imply that further participants are recruited.
Events or dropout can make observed follow-up end before planned
completion. `decision_time` is the look date for immediate success or
futility, and the analysis-ready date when a final analysis determines
success. Administrative data-cleaning delays are not modeled.

## Reusing trial paths for calibration

[`sim_horizon_paths()`](https://graemeleehickey.github.io/goldilocks/reference/sim_horizon_paths.md)
retains current-cohort predictive scores at every planned look and
supported horizon, maximum-cohort predictive scores at H_F, and final
analyses for every cohort and horizon.
[`apply_horizon_rules()`](https://graemeleehickey.github.io/goldilocks/reference/apply_horizon_rules.md)
selects the evidence for a candidate schedule and applies its
thresholds. Saving the individual scores, rather than only a predictive
probability at one `prob_ha`, permits recalibration of the
completed-data threshold as well.

Reuse can change the nondecreasing success-horizon schedule within the
retained looks and horizons, and `Fn`, `Sn`, `Qn`, and `prob_ha`. The
generating scenario, collection and final horizons, dropout model,
analysis method, priors, horizon-specific analysis settings, maximum
sample size, look sizes, and Monte Carlo counts remain fixed. Changes to
those assumptions require new paths.

Direct simulation with
[`sim_horizon_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_horizon_trials.md)
skips unreached interim calculations. With the same seed, settings, and
full set of `endpoint_horizons`, it agrees with complete-path simulation
followed by rule application. Explicit seeds preserve the caller’s
random-number state, and `trial_offset` supports consecutive batches.
Changing the retained set of horizons changes the allocation of random
numbers. Saved settings and calculation failures remain part of the path
object when it is written with
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html).

The ordinary futility measure is \Pr(\text{binding futility}) under the
simulated scenario. The counterfactual measure is

\Pr\\\text{binding futility and final success at }N\_{\max},H_F\\.

The latter is a joint probability over simulated trials, not a
conditional probability among futility stops and not the net loss of
power. It is reported as `futility_and_success_at_max`, separately from
`stop_futility`. A failed counterfactual final analysis makes the joint
outcome unavailable for a futility trial; available and missing
comparison counts accompany its summary. Failures needed to determine
the candidate’s actual adaptive outcome exclude that trial and are
reported separately.

Type I error must be assessed for the full adaptive design under the
relevant joint null scenarios across endpoints. Overall success can
occur at an earlier endpoint, so it does not necessarily constitute
rejection of the final-horizon null alone. Compare candidates using
paired trial-level outcomes from shared paths, and validate selected
settings using independent simulations. `N_trials` controls
operating-characteristic precision; `N_impute` and, where used, `N_mcmc`
also affect the implemented stopping procedure. [Calibrating stopping
rules from complete trial
paths](https://graemeleehickey.github.io/goldilocks/articles/calibrating-stopping-rules.md)
gives further guidance on comparisons and Monte Carlo uncertainty.
