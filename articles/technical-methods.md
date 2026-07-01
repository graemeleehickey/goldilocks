# Technical details of the Goldilocks design

## Vignette summary

The `goldilocks` package implements the Goldilocks adaptive sample-size
design of Broglio, Connor, and Berry (2014) for time-to-event endpoints.
This vignette outlines the technical details of the design, including
notation, the piecewise-exponential event-time model, Gamma posterior
updating, posterior **predictive probabilities**, interim decision
rules, final analysis options, and simulation-based calibration. The
package vignettes “Two-arm randomized trials”, “Bayesian
piecewise-exponential designs”, and “Single-arm designs with a
performance goal” provide more application-focused examples in R.

## 1. Design and notation

Consider a trial with maximum sample size $`N_{\max}`$, planned endpoint
time $`\tau`$, and interim sample-size selection analyses after

``` math
n_1 < n_2 < \cdots < n_L < N_{\max}
```

subjects have been enrolled. The package argument `N_total` corresponds
to $`N_{\max}`$, `end_of_study` to $`\tau`$, and `interim_look` to
$`(n_1,\ldots,n_L)`$.

Let $`Z_i \in \{0,1\}`$ denote treatment assignment for subject $`i`$,
matching the package data column `treatment`: $`Z_i = 1`$ for the
treatment arm and $`Z_i = 0`$ for the control arm. In a single-arm
design, all $`Z_i = 1`$ and there is no concurrent control. The package
assumes that randomization occurs at enrollment; in practice those times
can differ, but the distinction is not represented in the simulation
model. We therefore use enrollment time throughout. Let $`E_i`$ denote
enrollment time from first patient in, $`T_i^*`$ the true event time
from enrollment, $`C_i`$ the administrative or loss-to-follow-up
censoring time, and

``` math
T_i = \min(T_i^*, C_i), \qquad
  \delta_i = I(T_i^* \le C_i).
```

At the interim analysis held when $`n_\ell`$ subjects have been
enrolled, the analysis calendar time is $`E_{n_\ell}`$. Subject $`i`$’s
current observed follow-up is

``` math
u_{i\ell} = \max\{0, E_{n_\ell} - E_i\}.
```

For enrolled subjects, the interim data are masked to the information
available at this calendar time. If an event has not yet occurred by
$`u_{i\ell}`$, the subject is treated as censored at $`u_{i\ell}`$.
Subjects with $`i > n_\ell`$ are not yet enrolled and contribute only to
the maximum-sample-size prediction.

At each interim analysis the design estimates two predictive
probabilities:

1.  $`P_{n_\ell}`$, the probability of final success if accrual stops at
    the current enrolled sample size and all enrolled subjects complete
    follow-up;
2.  $`P_{\max,\ell}`$, the probability of final success if the trial
    continues to $`N_{\max}`$.

These are compared with thresholds $`S_\ell`$ and $`F_\ell`$,
corresponding to the package arguments `Sn` and `Fn`.

## 2. Event-time model

The time-to-event model is piecewise exponential. Let

``` math
0 = s_0 < s_1 < \cdots < s_{J-1}
```

be the vector of cut-points supplied through `cutpoints`. These define
intervals

``` math
[s_0, s_1), [s_1, s_2), \ldots, [s_{J-1}, \infty).
```

For treatment value $`z \in \{0,1\}`$, where $`z = 1`$ denotes the
treatment arm and $`z = 0`$ denotes the control arm, interval $`j`$ has
constant hazard $`\lambda_{zj}`$. The subject-level hazard is therefore

``` math
h_z(t) = \lambda_{zj}, \qquad s_{j-1} \le t < s_j,
```

where the last interval has no finite upper endpoint. With
`cutpoints = 0`, $`J = 1`$ and the model reduces to an ordinary
exponential model.

The cumulative hazard for treatment value $`z`$ at time $`t`$ is

``` math
H_z(t) =
  \sum_{j=1}^{J}
  \lambda_{zj}
  \{ \min(t, s_j) - s_{j-1} \}_+,
```

where $`s_J = \infty`$ and $`\{x\}_+ = \max(x,0)`$. The corresponding
survival and cumulative event probability are

``` math
S_z(t) = \exp\{-H_z(t)\}, \qquad
  p_z(t) = 1 - S_z(t).
```

The helper
[`prop_to_haz()`](https://graemeleehickey.github.io/goldilocks/reference/prop_to_haz.md)
solves the inverse problem used in simulation planning: given event
probabilities at one or more time points, it returns the piecewise
hazards that imply those probabilities. The helper
[`ppwe()`](https://graemeleehickey.github.io/goldilocks/reference/ppwe.md)
evaluates $`p_z(\tau)`$, and `haz_to_prop()` applies this transformation
to posterior hazard draws.

For observed follow-up $`(T_i, \delta_i, Z_i)`$, the
piecewise-exponential likelihood can be written in terms of
interval-specific event counts and exposure times. Define

``` math
d_{zj} = \sum_i I(Z_i = z)\delta_i I(T_i \in [s_{j-1}, s_j)),
```

and

``` math
y_{zj} =
  \sum_i I(Z_i = z)
  \{ \min(T_i, s_j) - s_{j-1} \}_+ I(T_i > s_{j-1}).
```

Up to factors not involving $`\lambda_{zj}`$, the likelihood
contribution for treatment value $`z`$ is

``` math
L_z(\boldsymbol{\lambda}_z; \mathcal{D})
  \propto
  \prod_{j=1}^{J}
  \lambda_{zj}^{d_{zj}} \exp(-\lambda_{zj}y_{zj}),
```

where $`\boldsymbol{\lambda}_z =
(\lambda_{z1},\ldots,\lambda_{zJ})^\top`$. This sufficient-statistic
form is what makes the Gamma posterior update available in closed form.

## 3. Posterior distribution of hazards

For each treatment value $`z`$ and interval $`j`$, `goldilocks` assumes
an independent Gamma prior

``` math
\lambda_{zj} \sim \operatorname{Gamma}(\alpha_0, \beta_0),
```

where $`\alpha_0`$ is the shape and $`\beta_0`$ is the rate. This
follows the [`stats::rgamma()`](https://rdrr.io/r/stats/GammaDist.html)
parameterization. The argument `prior = c(alpha0, beta0)` sets these two
hyperparameters, with default `prior = c(0.1, 0.1)`. The same prior is
applied to every treatment group and every piecewise interval. Separate
priors by treatment group or interval cannot currently be specified
through the package interface.

At an analysis, let $`d_{zj}`$ be the number of observed events for
treatment value $`z`$, interval $`j`$, and let $`y_{zj}`$ be the total
observed exposure time for that treatment value and interval.
Gamma-exponential conjugacy gives

``` math
\lambda_{zj} \mid \mathcal{D}
  \sim \operatorname{Gamma}(\alpha_0 + d_{zj}, \beta_0 + y_{zj}).
```

The package obtains $`(d_{zj}, y_{zj})`$ by splitting each subject’s
observed follow-up over the cut-point intervals. Posterior draws are
generated independently for each treatment group and interval. In early
interim analyses, later piecewise intervals may have no exposure. In
that case, the package propagates exposure and event counts from the
nearest non-empty interval within the same treatment group and emits a
warning. This prevents undefined posterior draws, but it should be
interpreted as a computational fallback rather than evidence about that
later interval.

The posterior density factorizes as

``` math
\pi(\boldsymbol{\lambda} \mid \mathcal{D})
  =
  \prod_z \prod_{j=1}^{J}
  \pi(\lambda_{zj} \mid d_{zj}, y_{zj}),
```

where each marginal factor is the Gamma distribution above. Posterior
predictive calculations integrate over this density rather than
conditioning on a single plug-in hazard estimate.

## 4. Predictive distribution for incomplete outcomes

At an interim analysis, subjects can be separated into three sets:

1.  enrolled subjects with complete endpoint information;
2.  enrolled subjects with partial follow-up; and
3.  future subjects not yet enrolled.

The first set contributes observed events and exposure to the posterior.
The second and third sets require prediction.

For an enrolled subject who is event-free through time $`u`$, a future
event time is drawn from the conditional piecewise-exponential
distribution

``` math
\Pr(T \le t \mid T > u)
  = \frac{F(t) - F(u)}{1 - F(u)}, \qquad t > u,
```

where $`F(t) = 1 - S(t)`$. Equivalently, if
$`U \sim \operatorname{Uniform}(0,1)`$, then

``` math
T = F^{-1}\{F(u) + U[1 - F(u)]\}.
```

This is implemented by
[`pwe_impute()`](https://graemeleehickey.github.io/goldilocks/reference/pwe_impute.md).
For future subjects in the maximum-sample-size calculation,
[`pwe_sim()`](https://graemeleehickey.github.io/goldilocks/reference/pwe_sim.md)
draws unconditional event times from the same piecewise-exponential
model. In the package implementation, these future event-time
imputations do not require explicitly simulating future enrollment times
at the interim look; the number of future subjects is determined by
$`N_{\max} - n_\ell`$.

Let $`\mathcal{D}_{\ell}^{\mathrm{obs}}`$ denote the data observed at
look $`\ell`$, and let $`\mathcal{D}^{\mathrm{mis}}`$ denote unobserved
event times and future censoring indicators. The posterior predictive
density is

``` math
p(\mathcal{D}^{\mathrm{mis}} \mid \mathcal{D}_{\ell}^{\mathrm{obs}})
  =
  \int
  p(\mathcal{D}^{\mathrm{mis}} \mid \boldsymbol{\lambda})
  \pi(\boldsymbol{\lambda} \mid \mathcal{D}_{\ell}^{\mathrm{obs}})
  d\boldsymbol{\lambda}.
```

This density is the mathematical target, not an object that the package
evaluates in closed form. In the implementation, `posterior()` draws
$`\boldsymbol{\lambda}`$, `impute_data()` draws completed outcomes
conditional on those hazards, and `test_stop_success()` repeats the
impute-and-analyze cycle to approximate the predictive probability by
Monte Carlo simulation.

## 5. Interim decision algorithm

At interim look $`\ell`$, the package first estimates the posterior of
the hazard parameters from the currently observable data. It then uses
Monte Carlo integration to approximate $`P_{n_\ell}`$ and
$`P_{\max,\ell}`$.

Let $`\psi(\mathcal{D})`$ be the final success indicator for a completed
analysis dataset:

``` math
\psi(\mathcal{D}) =
  I\{Q(\mathcal{D}) > c\},
```

where $`Q(\mathcal{D})`$ is the final analysis quantity and $`c`$ is the
success threshold. In package notation, $`c =`$`prob_ha`. The final
analysis quantity $`Q(\mathcal{D})`$ depends on the analysis method:

| Design setting | `method` | $`Q(\mathcal{D})`$ | Supported alternatives |
|----|----|----|----|
| Two-arm randomized trial | `logrank` | $`1-p(\mathcal{D})`$, where $`p(\mathcal{D})`$ is the traditional log-rank test P-value, with one-sided variants defined in Section 6.1 | `"less"`, `"greater"`, `"two.sided"` |
| Two-arm randomized trial | `cox` | $`1-p(\mathcal{D})`$, where $`p(\mathcal{D})`$ is the traditional Wald-test P-value, with one-sided variants defined in Section 6.1 | `"less"`, `"greater"`, `"two.sided"` |
| Two-arm randomized trial | `chisq` | $`1-p(\mathcal{D})`$, where $`p(\mathcal{D})`$ is the traditional chi-square test P-value | `"two.sided"` |
| Two-arm randomized trial | `bayes` | $`\Pr(\Delta < h_0 \mid \mathcal{D})`$ or $`\Pr(\Delta > h_0 \mid \mathcal{D})`$ | `"less"`, `"greater"` |
| Single-arm trial | `bayes` | $`\Pr(p_1(\tau) < h_0 \mid \mathcal{D})`$ or $`\Pr(p_1(\tau) > h_0 \mid \mathcal{D})`$ | `"less"`, `"greater"` |

For frequentist analyses, `prob_ha` is therefore a transformed P-value
threshold. For example, `prob_ha = 0.975` corresponds to a one-sided
$`\alpha = 0.025`$ rule. The value should be chosen during design
calibration to control the desired type I error rate across relevant
null scenarios.

### 5.1 Predictive probability at the current sample size

The current-sample-size predictive probability $`P_{n_\ell}`$ is
estimated by:

1.  draw one set of hazard parameters from the interim posterior;
2.  impute remaining endpoint outcomes for enrolled subjects only;
3.  apply the prespecified final analysis rule to the completed data;
    and
4.  record whether that completed trial is successful.

Formally,

``` math
P_{n_\ell}
  =
  \operatorname{E}\{
    \psi(\mathcal{D}_{n_\ell}^{\mathrm{comp}})
    \mid \mathcal{D}_{\ell}^{\mathrm{obs}}
  \},
```

where $`\mathcal{D}_{n_\ell}^{\mathrm{comp}}`$ is the completed dataset
formed from the $`n_\ell`$ enrolled subjects after imputing their
remaining follow-up. Equivalently,

``` math
P_{n_\ell}
  =
  \int
  \psi(\mathcal{D}_{n_\ell}^{\mathrm{obs}},
       \mathcal{D}_{n_\ell}^{\mathrm{mis}})
  p(\mathcal{D}_{n_\ell}^{\mathrm{mis}}
    \mid \mathcal{D}_{\ell}^{\mathrm{obs}})
  d\mathcal{D}_{n_\ell}^{\mathrm{mis}}.
```

Repeating the four-step procedure above for Monte Carlo replicate
$`m = 1,\ldots,M`$, where $`M =`$`N_impute`, gives

``` math
\widehat{P}_{n_\ell}
  = \frac{1}{M}\sum_{m=1}^{M} I\{\textrm{success in replicate } m\}.
```

If

``` math
\widehat{P}_{n_\ell} > S_\ell,
```

accrual is stopped for expected success. Enrolled subjects are still
followed to the planned final analysis time.

### 5.2 Predictive probability at the maximum sample size

The maximum-sample-size predictive probability $`P_{\max,\ell}`$ is
estimated similarly, except that the completed trial includes both
currently enrolled subjects and future subjects required to reach
$`N_{\max}`$. For each replicate, event times are imputed for the future
subjects, the completed dataset is analyzed, and success is recorded:

``` math
P_{\max,\ell}
  =
  \operatorname{E}\{
    \psi(\mathcal{D}_{N_{\max}}^{\mathrm{comp}})
    \mid \mathcal{D}_{\ell}^{\mathrm{obs}}
  \}.
```

``` math
\widehat{P}_{\max,\ell}
  = \frac{1}{M}\sum_{m=1}^{M}
    I\{\textrm{success at } N_{\max} \textrm{ in replicate } m\}.
```

If

``` math
\widehat{P}_{\max,\ell} < F_\ell,
```

the trial stops for futility. Otherwise, accrual continues to the next
interim look.

Thus the interim action at look $`\ell`$ can be represented as

``` math
A_\ell =
\begin{cases}
\textrm{stop accrual for expected success}, &
  \widehat{P}_{n_\ell} > S_\ell,\\
\textrm{stop for futility}, &
  \widehat{P}_{\max,\ell} < F_\ell,\\
\textrm{continue accrual}, &
  \textrm{otherwise}.
\end{cases}
```

The first look satisfying either stopping condition defines the adaptive
stopping look,

``` math
L^* = \inf\{\ell : \widehat{P}_{n_\ell} > S_\ell
  \textrm{ or } \widehat{P}_{\max,\ell} < F_\ell\},
```

with $`L^* = L + 1`$ if no interim stopping condition is met and the
design continues to $`N_{\max}`$.

## 6. Final analysis

The final analysis is conducted after accrual has stopped and the
relevant follow-up has completed for the enrolled cohort, subject to the
handling of loss to follow-up described below. The final rule supplies
the binary success indicator used inside the predictive probability
calculations.

### 6.1 Frequentist final tests

For `method = "logrank"`, success is based on a log-rank test. For
`method = "cox"`, success is based on the Wald test from a Cox
proportional hazards regression. For these methods, `goldilocks` stores
$`1-p`$ in `post_prob_ha`; this is not a posterior probability, but it
puts frequentist and Bayesian rules on a common “larger is stronger
evidence” scale. For example, a one-sided test at $`\alpha = 0.025`$
corresponds to `prob_ha = 0.975`.

For the log-rank option, let $`Z_{\mathrm{LR}}`$ denote the signed
log-rank statistic, with positive values corresponding to excess events
in the control arm under the package convention. Let $`p_{\mathrm{LR}}`$
denote the traditional two-sided log-rank test P-value. For the Cox
option, let $`\widehat{\eta}`$ be the estimated log hazard ratio for
treatment versus control and let

``` math
Z_{\mathrm{Cox}} = \frac{\widehat{\eta}}
  {\operatorname{se}(\widehat{\eta})}.
```

Here, lower treatment hazard corresponds to $`Z_{\mathrm{Cox}} < 0`$.
Let $`p_{\mathrm{Cox}}`$ denote the traditional two-sided Wald-test
P-value from the Cox model. The package uses the following
method-specific definitions of $`Q(\mathcal{D})`$:

| Method    | Alternative   | $`Q(\mathcal{D})`$             |
|-----------|---------------|--------------------------------|
| `logrank` | `"less"`      | $`\Phi(Z_{\mathrm{LR}})`$      |
| `logrank` | `"greater"`   | $`1 - \Phi(Z_{\mathrm{LR}})`$  |
| `logrank` | `"two.sided"` | $`1 - p_{\mathrm{LR}}`$        |
| `cox`     | `"less"`      | $`1 - \Phi(Z_{\mathrm{Cox}})`$ |
| `cox`     | `"greater"`   | $`\Phi(Z_{\mathrm{Cox}})`$     |
| `cox`     | `"two.sided"` | $`1 - p_{\mathrm{Cox}}`$       |
| `chisq`   | `"two.sided"` | $`1 - p_{\chi^2}`$             |

Here $`p_{\chi^2}`$ is the traditional chi-square test P-value. The
one-sided directions differ between the log-rank and Cox rows because of
the sign convention of the package’s log-rank statistic.

For `method = "chisq"`, the final event indicator is compared between
arms using a chi-square test. This discards event-time information and
is available only for two-arm designs with `alternative = "two.sided"`.

### 6.2 Bayesian final test

For `method = "bayes"`, posterior hazard draws are mapped to cumulative
event probabilities at $`\tau`$. In a two-arm design the treatment
effect is

``` math
\Delta = p_1(\tau) - p_0(\tau),
```

where $`p_1(\tau)`$ is the treatment-arm event probability and
$`p_0(\tau)`$ is the control-arm event probability. The effect is on the
event scale, not the survival scale. For an adverse event, benefit
usually means $`\Delta < 0`$.

Because $`p_a(\tau) = 1 - \exp\{-H_a(\tau)\}`$, posterior draws of
$`\Delta`$ are obtained by transforming posterior draws of
$`\boldsymbol{\lambda}_1`$ and $`\boldsymbol{\lambda}_0`$:

``` math
\Delta^{(b)}
  =
  \left[1 - \exp\{-H_1^{(b)}(\tau)\}\right]
  -
  \left[1 - \exp\{-H_0^{(b)}(\tau)\}\right],
  \qquad b = 1,\ldots,B.
```

The Monte Carlo estimate of the posterior probability for
`alternative = "less"` is

``` math
\widehat{\Pr}(\Delta < h_0 \mid \mathcal{D})
  =
  \frac{1}{B}\sum_{b=1}^{B} I(\Delta^{(b)} < h_0),
```

where $`B =`$`N_mcmc`.

With `alternative = "less"`, success is declared when

``` math
\Pr(\Delta < h_0 \mid \mathcal{D}) > \texttt{prob_ha}.
```

With `alternative = "greater"`, success is declared when

``` math
\Pr(\Delta > h_0 \mid \mathcal{D}) > \texttt{prob_ha}.
```

The Bayesian final test is one-sided in the package;
`alternative = "two.sided"` is not supported.

In a single-arm design there is no $`p_0(\tau)`$. The estimand becomes
$`p_1(\tau)`$, and $`h_0`$ is an external benchmark event probability.
In clinical-trial terminology this benchmark is often called a
performance goal (PG) or objective performance criterion (OPC).
Consequently, single-arm designs in `goldilocks` require
`method = "bayes"`.

### 6.3 Loss to follow-up at the final analysis

Interim predictions impute outcomes that are not yet known. At the final
analysis, `imputed_final` controls whether subjects lost to follow-up
are also imputed.

If `imputed_final = TRUE`, the final analysis mirrors the interim
predictive framework: each imputed completed dataset is analyzed,
yielding a value of the success probability scale $`Q(\mathcal{D})`$
(posterior probability for `method = "bayes"` or transformed P-value for
frequentist methods), and `post_prob_ha` is the average of those values
across imputations. If `imputed_final = FALSE`, the final analysis uses
observed right-censored data for methods that can handle censoring
(`logrank`, `cox`, and `bayes`). For `chisq`, lost-to-follow-up subjects
are excluded when `imputed_final = FALSE` because the chi-square test
has no mechanism for right-censored observations.

The loss-to-follow-up mechanism in the simulator is non-informative.
Designs where dropout may depend on prognosis should be assessed with
sensitivity analyses outside the default data-generating mechanism.

## 7. Operating characteristics

A Goldilocks design is calibrated by simulation. A single simulated
trial describes one possible path; the design is characterized by
repeated simulation over clinically relevant scenarios.

For a candidate design,
[`sim_trials()`](https://graemeleehickey.github.io/goldilocks/reference/sim_trials.md)
repeatedly calls
[`survival_adapt()`](https://graemeleehickey.github.io/goldilocks/reference/survival_adapt.md)
and
[`summarise_sims()`](https://graemeleehickey.github.io/goldilocks/reference/summarise_sims.md)
estimates:

- power, or type I error under a null scenario;
- probability of stopping for expected success;
- probability of stopping for futility;
- probability of reaching $`N_{\max}`$;
- mean and standard deviation of enrolled sample size; and
- probability of stopping for expected success but failing at the final
  analysis.

Let $`R = 1,\ldots,R_{\max}`$ index simulated trials under a scenario
$`\theta`$, where $`\theta`$ denotes the data-generating parameters such
as control hazard, treatment hazard, accrual rate, loss-to-follow-up
rate, and follow-up duration. The trial-level random variables are:

| Symbol  | Meaning                                                         |
|---------|-----------------------------------------------------------------|
| $`N_R`$ | enrolled sample size in simulated trial $`R`$                   |
| $`E_R`$ | indicator that trial $`R`$ stopped accrual for expected success |
| $`F_R`$ | indicator that trial $`R`$ stopped for futility                 |
| $`Z_R`$ | final success indicator in trial $`R`$                          |

Let $`\Theta_0`$ denote the null parameter space, i.e. the set of
data-generating scenarios in which the treatment does not satisfy the
alternative hypothesis. For a two-arm superiority trial this includes
scenarios with no beneficial treatment effect; in practice, it should be
explored across plausible nuisance parameters such as control event
rates and accrual rates. The main operating characteristics are

``` math
\operatorname{Power}(\theta)
  = \Pr_\theta(Z_R = 1, F_R = 0),
```

``` math
\operatorname{Type\ I\ error}(\theta_0)
  = \Pr_{\theta_0}(Z_R = 1, F_R = 0), \qquad \theta_0 \in \Theta_0,
```

``` math
\Pr_\theta(\textrm{stop for expected success}) = \Pr_\theta(E_R = 1),
```

``` math
\Pr_\theta(\textrm{stop for futility}) = \Pr_\theta(F_R = 1),
```

and

``` math
\operatorname{E}_\theta(N_R), \qquad
  \operatorname{Var}_\theta(N_R).
```

The `stop_and_fail` summary estimates

``` math
\Pr_\theta(E_R = 1, Z_R = 0),
```

which is the probability that accrual stops for expected success but the
final analysis does not meet the success criterion.

The basic workflow is:

``` r

out <- sim_trials(
  hazard_treatment = ht,
  hazard_control   = hc,
  cutpoints        = cutpoints,
  N_total          = N_total,
  lambda           = lambda,
  lambda_time      = lambda_time,
  interim_look     = interim_look,
  end_of_study     = end_of_study,
  prior            = prior,
  Fn               = Fn,
  Sn               = Sn,
  prob_ha          = prob_ha,
  N_impute         = N_impute,
  N_mcmc           = N_mcmc,
  N_trials         = N_trials,
  method           = method)

summarise_sims(out$sims)
```

Broglio et al. emphasize that type I error for this class of adaptive
design should be examined across the null space, not only at one
convenient null scenario. For time-to-event endpoints, the relevant null
space includes plausible control event rates and accrual rates. Accrual
rate is especially important because rapid enrollment can leave little
endpoint information available at interim looks, increasing the
uncertainty in both $`\widehat{P}_{n_\ell}`$ and
$`\widehat{P}_{\max,\ell}`$. The follow-up period after accrual stops
also affects operating characteristics because it determines how much
additional information is observed before the final analysis.

The Monte Carlo sizes `N_impute`, `N_mcmc`, and `N_trials` should be
chosen so that simulation error is small relative to the decision being
made. Small values are appropriate for examples and package tests, but
final design calibration usually requires many more trials and
imputations.

## 8. Threshold selection

The thresholds $`S_\ell`$ and $`F_\ell`$ may be constant across looks or
may vary by look. They interact with the final analysis threshold,
number and timing of looks, endpoint delay, accrual rate, loss to
follow-up, prior distribution, and maximum sample size.

If type I error is too high, possible remedies include increasing
`prob_ha`, increasing the expected-success thresholds, reducing the
number of looks, or altering follow-up requirements. If power is too
low, the maximum sample size, futility threshold, expected-success
threshold, or final analysis threshold may need reconsideration. Each
change should be rechecked under null and alternative scenarios.

The quantity `stop_and_fail` is particularly useful when tuning
$`S_\ell`$. It estimates how often a trial stops accrual for expected
success but does not meet the final success criterion after follow-up is
complete. If this is too large, the expected-success threshold is
usually too permissive for the amount of uncertainty present at interim
looks.

## 9. Relation to group-sequential designs

A group-sequential design usually indexes interim analyses by
information, such as the number of observed events, and may stop
immediately for efficacy when a boundary is crossed. A Goldilocks design
indexes sample-size selection analyses by enrolled sample size and
explicitly incorporates future follow-up of the currently enrolled
cohort.

This distinction matters for delayed outcomes. A trial may enroll many
subjects before accumulating enough events for an event-driven interim
analysis. Goldilocks uses the partial information available during
accrual to decide whether additional subjects are needed, while
preserving a single preplanned final analysis after the enrolled
subjects have completed follow-up.

## 10. Package-specific scope

The time-to-event example in Broglio et al. uses a Gamma-exponential
prediction model. `goldilocks` extends this to a piecewise-exponential
model, so the hazard may change at prespecified cut-points. This can be
useful when there is a clinically plausible early-risk period, but each
additional interval adds parameters and can make interim posteriors more
diffuse.

The package also supports single-arm Bayesian Goldilocks designs by
replacing the concurrent control with an external benchmark $`h_0`$,
often referred to as a performance goal (PG) or objective performance
criterion (OPC). This is convenient for early-phase, rare-disease, or
proof-of-concept settings, but validity then depends on the benchmark
being transportable to the enrolled population.

## References

Broglio KR, Connor JT, Berry SM. Not too big, not too small: a
Goldilocks approach to sample size selection. *Journal of
Biopharmaceutical Statistics*, 2014; **24(3)**: 685-705.
<doi:10.1080/10543406.2014.888569>.

U.S. Food and Drug Administration. *Adaptive Design Clinical Trials for
Drugs and Biologics Guidance for Industry*. December 2019.
<https://www.fda.gov/regulatory-information/search-fda-guidance-documents/adaptive-design-clinical-trials-drugs-and-biologics-guidance-industry>.

U.S. Food and Drug Administration. *Guidance for the Use of Bayesian
Statistics in Medical Device Clinical Trials*. February 5, 2010.
<https://www.fda.gov/media/71512/download>.
