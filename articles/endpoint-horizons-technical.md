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
prediction. All participants have planned outcome collection through H_F
from enrollment, including participants recruited before the endpoint
changes. `final_horizon` therefore specifies both the eventual endpoint
and collection duration; no separate collection-horizon argument is
needed.

| Quantity | R argument | Interpretation |
|:---|:---|:---|
| n\_\ell, H\_\ell | `endpoint_schedule` | Data frame with columns `N` and `horizon`; one row per interim |
| N\_{\max} | `N_total` | Maximum enrollment; not an additional interim look |
| H_F | `final_horizon` | Maximum-N endpoint, futility prediction horizon, and planned collection duration |
| Dropout reference time | `dropout_horizon` | Time at which `prop_loss` specifies the dropout probability |

The schedule is applied before evaluating each reached look. It does not
use an interim estimate or predictive probability to choose whether to
extend the endpoint. A terminal stopping decision prevents later
scheduled changes from taking effect. If the trial continues to maximum
enrollment, its final analysis uses H_F after endpoint ascertainment;
reaching N\_{\max} does not itself trigger an interim analysis.

> **An expected-success stop fixes the cohort and its endpoint.** If
> accrual stops at N=300 using a 12-month endpoint, the subsequent final
> analysis uses those 300 participants at 12 months. If that analysis
> fails, the trial remains unsuccessful; it does not resume enrollment
> or switch to an 18-month analysis to establish success. Planned
> collection through 18 months does not change this rule. Expected
> success is a reason to stop accrual, not a declaration of final
> success.

## Two data cuts at the same calendar date

Let E_i be enrollment time, T_i the event time from enrollment, and D_i
the independent dropout time. The calendar date of look \ell is t\_\ell
= E\_{n\_\ell}. For an enrolled participant and an analysis horizon H,
the available data are

X\_{i\ell}(H) = \min\\T_i, D_i, H, t\_\ell-E_i\\, \qquad
\delta\_{i\ell}(H) = I\\T_i \le \min(D_i,H,t\_\ell-E_i)\\.

The success data set uses H=H\_\ell; the futility data set uses H=H_F.
Here, a **data cut** means the observations retained after restricting
both calendar availability and each participant’s follow-up duration.
Participants not yet enrolled contribute no observed events or exposure
to either predictive model. Their predictive outcomes are generated
through H_F, only for the futility assessment of success at maximum N.
The implementation gives a participant enrolled exactly at the look a
negligible positive exposure to accommodate the survival calculations.

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

**If t_1 \< H \< t_2, distinguish the hazard parameter from the
observations used to estimate it.** The master model still specifies one
constant hazard on (t_1,t_2\]. For success prediction at H, its
likelihood uses only events and exposure in (t_1,H\] that are available
at the interim date. The original Gamma prior for that hazard is
retained without rescaling. Thus, with prior \operatorname{Gamma}(a,b),
the update is \operatorname{Gamma}\\a+d(t_1,H),b+y(t_1,H)\\. The fitted
interval is displayed as (t_1,H\] because the success prediction ends
there; this does not introduce a new underlying hazard change at H.

For example, with cutpoints at months 6 and 15 and success endpoint 12,
the second hazard is updated using observations from months 6 to 12. An
observed event at month 14 supplies no event or post-12-month exposure
to that success update. The separate futility model targeting month 18
can use the event and exposure through month 14 in its (6,15\] interval,
provided they are already observed on the look date.

An alternative would estimate the same hazard using all available
observations in (t_1,t_2\], then predict only through H. This borrows
information beyond the success endpoint under the constant-hazard
assumption and may improve precision if that assumption is appropriate.
It also allows post-H observations to change success stopping. **The
current functions deliberately use the first approach: no observations
beyond H inform success prediction, even within a shared hazard
interval.** Futility uses its own observations through H_F.

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
imputation at the actual final analysis when the selected method and
empty-interval policy require it.

The optional `horizon_settings` data frame specifies horizon-specific
`h0`, `rmst_tau`, or `prior_bin`. Without overrides, the common null
value and beta prior are used at every horizon, while the RMST
restriction equals the endpoint horizon. This allows, for example,
different single-arm performance goals at different endpoint times. The
analysis method and alternative direction remain common across horizons.
The [applied
vignette](https://graemeleehickey.github.io/goldilocks/articles/endpoint-horizons-example.html#optional-endpoint-specific-margins)
shows how to specify different binary-risk margins at 12 and 18 months.

When success and futility use different endpoint times, they fit
different observed data sets and make separate posterior and imputation
draws. For example, the first completed data set for 12-month success
and the first completed data set for 18-month futility are two
predictive simulations, not the same simulated future followed for
different durations. Their draw numbers are bookkeeping labels, not
matched patient histories. The two sets estimate separate marginal
predictive probabilities; they cannot be paired draw by draw to estimate
a joint event. The underlying trial history used to assess the eventual
trial outcomes is generated once and is shared across horizons.

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

The endpoint remains fixed after an expected-success stop, as
highlighted above. The futility probability concerns the final analysis
at N\_{\max} and H_F; it is not the probability of declaring success at
any remaining interim under the entire adaptive policy.

Final analysis with dropout follows the method-specific choices for
`imputed_final`. For `method = "bayes-surv"` with
`empty_interval = "prior"`, both flags use the observed-data posterior
at the selected endpoint, retaining observed events and censored
exposure. Averaging completed-data posterior probabilities and mean
effects under the same model and prior targets these summaries directly,
so final imputation is redundant; `N_mcmc` controls final precision,
while `N_impute` still controls interim prediction. For other settings,
`FALSE` retains right-censoring in survival analyses and excludes
incomplete endpoint statuses in binary analyses; with `TRUE`, supported
methods impute missing endpoint outcomes. Independent dropout does not
by itself make complete-case binary analysis unbiased. These choices and
the restrictions on log-rank and Farrington-Manning final imputation are
described in the [technical methods
vignette](https://graemeleehickey.github.io/goldilocks/articles/technical-methods.html#loss-to-follow-up-at-the-final-analysis).

For the final application of this conjugate Bayesian survival analysis
at the selected endpoint, use substantially more than the default 1,000
draws: 50,000 is a starting point, with a larger prespecified budget,
such as 500,000, near the success threshold. Assess Monte Carlo
uncertainty relative to that threshold. `N_mcmc` also controls posterior
draws within interim predictions at each endpoint, so increasing it
increases interim computation; there is no separate final-only draw
argument. Assess operating characteristics with the intended settings,
including the selected endpoint’s final draw count.

## Shared histories, dropout, and timing

Event times and dropout times are generated once per simulated trial
through `final_horizon`. Final cohort data are re-censored at each
endpoint; changing the endpoint does not redraw either process.
`generation_cutpoints` and the arm-specific generating hazards must
describe the entire collection period, and need not match the analysis
partition.

If the arm-specific dropout probability is q_z by reference time H_D,
then

D_i\mid Z_i=z \sim \operatorname{Exponential}(\eta_z),\qquad
\eta_z=-\log(1-q_z)/H_D.

Here, H_D= `dropout_horizon`. This hazard remains fixed across endpoint
times. With collection through H_F, the probability of dropout by that
time is 1-\exp(-\eta_z H_F); the proportion actually censored by dropout
can be smaller because events may occur first.

### Analysis and decision dates

`analysis_ready_time` is the calendar date by which the event and
censoring observations used in the selected endpoint’s final analysis
are available. `decision_time` is when the trial’s success or failure is
determined. Both are measured from first enrollment and exclude
data-cleaning delays.

| Stopping rule | `decision_time` |
|:---|:---|
| Expected success | The subsequent final analysis date, `analysis_ready_time` |
| Immediate success or futility | The interim look date |
| Maximum enrollment | The final analysis date, `analysis_ready_time` |

For example, suppose recruitment stops for expected success at calendar
month 30 using a 12-month endpoint. If the last participant remains
event-free and observed for 12 months, both the final analysis and
decision occur at month 42. A futility stop at month 30 instead
determines failure immediately. Continued follow-up does not change that
decision.

## Reusing trial paths for calibration

### How the schedule determines saved horizons

Both simulation functions take `endpoint_schedule` and `final_horizon`.
Look sizes come from `endpoint_schedule$N`; the endpoint times to
calculate are derived as
`sort(unique(c(endpoint_schedule$horizon, final_horizon)))`. The final
horizon is included even if every interim uses a shorter endpoint. There
is no separate `interim_look` or `endpoint_horizons` argument in these
functions.

For the schedule (100,12), (200,12), (300,12), (400,18),
[`sim_horizon_paths()`](https://graemeleehickey.github.io/goldilocks/reference/sim_horizon_paths.md)
evaluates current-cohort success at both 12 and 18 months at each of the
four looks. It also saves maximum-cohort futility predictions at H_F=18
and final analyses for each cohort and endpoint.
[`apply_horizon_rules()`](https://graemeleehickey.github.io/goldilocks/reference/apply_horizon_rules.md)
then selects one saved success endpoint per look. This permits
comparison of a switch at N=300 with a switch at N=400 without new
simulations. A candidate requiring an endpoint time absent from the
original schedule and final horizon needs new paths.

Saving the individual completed-data scores, rather than only a
predictive probability at one `prob_ha`, also permits recalibration of
that success threshold. Reuse can change the nondecreasing
success-horizon schedule within the retained looks and horizons, and
`Fn`, `Sn`, `Qn`, and `prob_ha`. The generating scenario, final horizon,
dropout model, analysis method, priors, horizon-specific analysis
settings, maximum sample size, look sizes, and Monte Carlo counts remain
fixed.

Direct simulation with
[`sim_horizon_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_horizon_trials.md)
skips unreached interim calculations. With the same seed, settings, and
set of derived horizons, it agrees with complete-path simulation
followed by rule application. Changing only the switch point, while
retaining the same endpoint times, preserves this agreement. Removing an
endpoint time entirely changes the allocation of random numbers, so
direct simulation and replay from the larger horizon set need not match
trial by trial. Saved settings and calculation failures remain part of
the path object when written with
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html); previously saved
paths can still be evaluated with
[`apply_horizon_rules()`](https://graemeleehickey.github.io/goldilocks/reference/apply_horizon_rules.md).

## Interpreting the two futility measures

Both measures are also available for fixed-horizon designs through
[`sim_trial_paths()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trial_paths.md)
and
[`apply_stopping_rules()`](https://graemeleehickey.github.io/goldilocks/reference/apply_stopping_rules.md);
they are not specific to endpoint extensions. See [Calibrating stopping
rules from complete trial
paths](https://graemeleehickey.github.io/goldilocks/articles/calibrating-stopping-rules.html#two-distinct-futility-measures).

Suppose 1,000 evaluable trials are simulated under an alternative, 200
stop for futility, and 40 of those 200 would have succeeded if forced to
continue to N\_{\max} and H_F. The two reported measures answer
different questions:

| Measure | Question | Example |
|:---|:---|---:|
| `stop_futility` | What fraction of all trials stop for futility? | 200/1000=20\\ |
| `futility_and_success_at_max` | What fraction of all trials stop for futility but would succeed at maximum N and H_F? | 40/1000=4\\ |

The second measure uses the same simulated trial history, continued
despite its futility stop. Its mathematical definition is

\Pr\\\text{binding futility and final success at }N\_{\max},H_F\\.

Among futility stops alone, the corresponding fraction would be
40/200=20\\. That conditional fraction is not what
`futility_and_success_at_max` reports.

**When futility is the only early stopping rule, the 4% is exactly a
4-percentage-point loss of power.** The comparison is with no early
stopping (`Fn = 0`, `Sn = Qn = 1`): all 1,000 trials reach the maximum-N
final analysis. Adding binding futility classifies the 40 trials that
would otherwise succeed as unsuccessful. This assumes matched trial
histories, analysis settings, and complete calculations. When other
stopping rules are present, such as stopping for expected success, the
interpretation is more complex.

## Type I error across endpoint times

The design can declare success at 12 months or at 18 months. A starting
null scenario is therefore one in which the treatment has no benefit at
either endpoint. For example, identical control and treatment hazards
throughout follow-up give identical event risks at both times. Simulate
the entire procedure, including the endpoint switches and all stopping
rules, and count any declaration of success. Under this scenario, that
fraction is the overall type I error. Calibrating each endpoint
separately does not establish control for the procedure that can stop at
either one.

The null assumptions must match the hypotheses used by the chosen
analysis method. For binary risk comparisons, equal risks at 18 months
alone do not imply equal risks at 12 months: an early treatment benefit
can disappear later. A success declared at 12 months then concerns the
12-month hypothesis; it does not establish benefit at 18 months. If the
design is intended to support endpoint-specific claims with strong
familywise error control, assess plausible configurations where only
some endpoint null hypotheses are true as well. In those configurations,
count success at a truly null endpoint, not every success; the overall
`power` column also includes legitimate successes at endpoints with
benefit. The [applied
example](https://graemeleehickey.github.io/goldilocks/articles/endpoint-horizons-example.html#benefit-at-12-months-but-no-benefit-at-18-months)
uses 12-month risks of 22% versus 30%, but an identical 18-month risk of
40.5%, and reports successes separately by endpoint. The distinction
between errors at individual endpoints and overall error is discussed in
FDA (2022, Sections II.B-C).

For a Cox analysis under a correctly specified proportional-hazards
model, with the same null margin at every look, the treatment hazard
ratio is common across follow-up horizons. Extending from 12 to 18
months then changes the available information, rather than the
underlying hazard ratio. In this setting, the differing endpoint nulls
illustrated by the binary-risk example are less of a concern. Multiple
opportunities to stop for success still require assessment of the full
stopping procedure.

Without proportional hazards, this simplification does not apply. A
single Cox hazard ratio summarizes effects over the included follow-up,
so extending the horizon can change the effect being summarized (Hernán,
2010). The log-rank test compares event-time distributions over the
included period, not just the event risks at its endpoint. Equal risks
at month 18 therefore do not establish its null hypothesis when the
earlier survival curves differ; early benefit and later harm can also
offset one another in the test. Switching to log-rank or Cox analysis
does not, by itself, remove the need to specify the relevant null
scenario.

## References

U.S. Food and Drug Administration. *Multiple Endpoints in Clinical
Trials: Guidance for Industry*. October 2022. [Guidance
document](https://www.fda.gov/media/162416/download).

Hernán MA. The hazards of hazard ratios. *Epidemiology*.
2010;21(1):13-15.
[doi:10.1097/EDE.0b013e3181c1ea43](https://doi.org/10.1097/EDE.0b013e3181c1ea43).
