# Technical details of the Goldilocks design

## Vignette summary

The `goldilocks` package implements the Goldilocks adaptive sample-size
design of Broglio, Connor, and Berry (2014) for time-to-event endpoints.
This vignette outlines the technical details of the design, including
notation, the continuous-time enrollment process, the
piecewise-exponential event-time model, Gamma posterior updating,
posterior **predictive probabilities**, interim decision rules, final
analysis options, and simulation-based calibration. The package
vignettes “Two-arm randomized trials”, “Bayesian piecewise-exponential
designs”, and “Single-arm designs with a performance goal” provide more
application-focused examples in R.

## 1. Design and notation

Consider a trial with maximum sample size N\_{\max}, planned endpoint
time \tau, and interim sample-size selection analyses after

n_1 \< n_2 \< \cdots \< n_L \< N\_{\max}

subjects have been enrolled. The package argument `N_total` corresponds
to N\_{\max}, `end_of_study` to \tau, and `interim_look` to
(n_1,\ldots,n_L).

Let Z_i \in \\0,1\\ denote treatment assignment for subject i, matching
the package data column `treatment`: Z_i = 1 for the treatment arm and
Z_i = 0 for the control arm. In a single-arm design, all Z_i = 1 and
there is no concurrent control. The package assumes that randomization
occurs at enrollment; in practice those times can differ, but the
distinction is not represented in the simulation model. We therefore use
enrollment time throughout. Let E_i denote enrollment time from first
patient in, T_i^\* the true event time from enrollment, C_i the
administrative or loss-to-follow-up censoring time, and

T_i = \min(T_i^\*, C_i), \qquad \delta_i = I(T_i^\* \le C_i).

At the interim analysis held when n\_\ell subjects have been enrolled,
the analysis calendar time is E\_{n\_\ell}. Subject i’s current observed
follow-up is

u\_{i\ell} = \max\\0, E\_{n\_\ell} - E_i\\.

For enrolled subjects, the interim data are masked to the information
available at this calendar time. If an event has not yet occurred by
u\_{i\ell}, the subject is treated as censored at u\_{i\ell}. Subjects
with i \> n\_\ell are not yet enrolled and contribute only to the
maximum-sample-size prediction.

At each interim analysis the design estimates two predictive
probabilities:

1.  P\_{n\_\ell}, the probability of final success if accrual stops at
    the current enrolled sample size and all enrolled subjects complete
    follow-up;
2.  P\_{\max,\ell}, the probability of final success if the trial
    continues to N\_{\max}.

These are compared with thresholds S\_\ell and F\_\ell, corresponding to
the package arguments `Sn` and `Fn`.

## 2. Continuous-time enrollment process

The simulation model places the first enrolled subject at time zero.
This is a first-patient-in origin, not an earlier protocol-approval,
site-activation, or recruitment-opening date. Let E_1 = 0 and let E_i
for i \> 1 denote the calendar time of the ith enrollment measured from
first patient in.

Enrollment after the first patient follows a non-homogeneous Poisson
process with a piecewise-constant intensity. Write the K internal
enrollment-rate knots as

0 \< a_1 \< a_2 \< \cdots \< a_K,

and define a_0 = 0 and a\_{K+1} = \infty. These internal knots are
supplied as `lambda_time`; zero is implicit and is not included. There
are K+1 positive rates

\boldsymbol{\rho} = (\rho_1,\ldots,\rho\_{K+1}),

supplied as `lambda`, with

\rho(t) = \rho_j, \qquad a\_{j-1} \le t \< a_j.

Consequently, `length(lambda)` must be exactly
`length(lambda_time) + 1`. A constant enrollment rate is represented by
`lambda_time = NULL` and a scalar `lambda`. The final rate continues
beyond the last knot until the requested `N_total` is reached;
[`enrollment()`](https://graemeleehickey.github.io/goldilocks/reference/enrollment.md)
does not impose a finite recruitment horizon.

The cumulative enrollment intensity is

A(t) = \int_0^t \rho(u)\\du = \sum\_{j=1}^{K+1}
\rho_j\\\min(t,a_j)-a\_{j-1}\\\_+.

To generate the process exactly, draw independent variables
X_2,\ldots,X\_{N\_{\max}} \sim \operatorname{Exponential}(1) and form

Q_i = \sum\_{k=2}^{i} X_k.

The enrollment times are then obtained through the inverse cumulative
intensity,

E_1 = 0, \qquad E_i = A^{-1}(Q_i), \quad i=2,\ldots,N\_{\max}.

More explicitly, if A(a\_{j-1}) \le Q_i \< A(a_j), then

E_i = a\_{j-1} + \frac{Q_i-A(a\_{j-1})}{\rho_j}.

This time-rescaling construction is exact: after the patient anchored at
zero, the number of arrivals in any interval (s,t\] is Poisson with mean
A(t)-A(s), and counts over disjoint intervals are independent. Under a
constant rate \rho, the successive gaps are independent
\operatorname{Exponential}(\rho) variables and E_n \sim
\operatorname{Gamma}(n-1,\text{rate}=\rho). In particular,
\operatorname{E}(E_n)=(n-1)/\rho. The fixed first patient accounts for
the n-1, rather than n, random gaps.

This differs from generating Poisson counts in unit-time bins and adding
uniform jitter afterward. Binning is sensitive to the arbitrary width of
a time unit and cannot represent a rate change inside a bin.
Cumulative-intensity inversion handles integer and fractional knots
identically and requires work proportional to the requested number of
subjects rather than to the elapsed number of empty bins.

For example,

``` r

enrollment(
  lambda = c(2, 5, 8),
  lambda_time = c(3.5, 9),
  N_total = 100
)
```

uses rate 2 over \[0,3.5), rate 5 over \[3.5,9), and rate 8 thereafter.
All rates are enrollments per common time unit. `lambda_time`,
enrollment times, event times, hazard `cutpoints`, and `end_of_study`
should therefore all use the same unit, such as days or months.

Although `lambda_time` and hazard `cutpoints` share the same
internal-knot API, they operate on different clocks. Enrollment knots
are trial-calendar times from first patient in. Hazard cutpoints are
follow-up times from each individual subject’s enrollment. Their values
and lengths are unrelated unless the scientific design happens to make
them coincide.

The package does not currently model site activations, site-specific
random rates, pauses, recruitment caps, or uncertainty in the supplied
rates. Those operational features require a richer site-level accrual
model; `lambda` represents the trial-level rate schedule assumed for a
simulation scenario.

## 3. Event-time model

The time-to-event model is piecewise exponential. Let

0 = s_0 \< s_1 \< \cdots \< s\_{J-1}.

The interior cutpoints (s_1,\ldots,s\_{J-1}) are supplied through
`cutpoints`; the initial boundary s_0 = 0 is implicit. These define
intervals

\[s_0, s_1), \[s_1, s_2), \ldots, \[s\_{J-1}, \infty).

For treatment value z \in \\0,1\\, where z = 1 denotes the treatment arm
and z = 0 denotes the control arm, interval j has constant hazard
\lambda\_{zj}. The subject-level hazard is therefore

h_z(t) = \lambda\_{zj}, \qquad s\_{j-1} \le t \< s_j,

where the last interval has no finite upper endpoint. With
`cutpoints = NULL`, J = 1 and the model reduces to an ordinary
exponential model.

The cumulative hazard for treatment value z at time t is

H_z(t) = \sum\_{j=1}^{J} \lambda\_{zj} \\ \min(t, s_j) - s\_{j-1} \\\_+,

where s_J = \infty and \\x\\\_+ = \max(x,0). The corresponding survival
and cumulative event probability are

S_z(t) = \exp\\-H_z(t)\\, \qquad p_z(t) = 1 - S_z(t).

The helper
[`prop_to_haz()`](https://graemeleehickey.github.io/goldilocks/reference/prop_to_haz.md)
solves the inverse problem used in simulation planning: given event
probabilities at one or more time points, it returns the piecewise
hazards that imply those probabilities. The helper
[`ppwe()`](https://graemeleehickey.github.io/goldilocks/reference/ppwe.md)
evaluates p_z(\tau), and `haz_to_prop()` applies this transformation to
posterior hazard draws.

For observed follow-up (T_i, \delta_i, Z_i), the piecewise-exponential
likelihood can be written in terms of interval-specific event counts and
exposure times. Define

d\_{zj} = \sum_i I(Z_i = z)\delta_i I(T_i \in \[s\_{j-1}, s_j)),

and

y\_{zj} = \sum_i I(Z_i = z) \\ \min(T_i, s_j) - s\_{j-1} \\\_+ I(T_i \>
s\_{j-1}).

Up to factors not involving \lambda\_{zj}, the likelihood contribution
for treatment value z is

L_z(\boldsymbol{\lambda}\_z; \mathcal{D}) \propto \prod\_{j=1}^{J}
\lambda\_{zj}^{d\_{zj}} \exp(-\lambda\_{zj}y\_{zj}),

where \boldsymbol{\lambda}\_z =
(\lambda\_{z1},\ldots,\lambda\_{zJ})^\top. This sufficient-statistic
form is what makes the Gamma posterior update available in closed form.

## 4. Posterior distribution of hazards

For each treatment value z and interval j, `goldilocks` assumes an
independent Gamma prior

\lambda\_{zj} \sim \operatorname{Gamma}(\alpha_0, \beta_0),

where \alpha_0 is the shape and \beta_0 is the rate. This follows the
[`stats::rgamma()`](https://rdrr.io/r/stats/GammaDist.html)
parameterization. The argument `prior = c(alpha0, beta0)` sets these two
hyperparameters, with default `prior = c(0.1, 0.1)`. The same prior is
applied to every treatment group and every piecewise interval. Separate
priors by treatment group or interval cannot currently be specified
through the package interface.

At an analysis, let d\_{zj} be the number of observed events for
treatment value z, interval j, and let y\_{zj} be the total observed
exposure time for that treatment value and interval. Gamma-exponential
conjugacy gives

\lambda\_{zj} \mid \mathcal{D} \sim \operatorname{Gamma}(\alpha_0 +
d\_{zj}, \beta_0 + y\_{zj}).

The package obtains (d\_{zj}, y\_{zj}) by splitting each subject’s
observed follow-up over the cut-point intervals. Posterior draws are
generated independently for each treatment group and interval. In early
interim analyses, later piecewise intervals may have no exposure. The
`empty_interval` argument controls the policy for these intervals:

- `empty_interval = "propagate"` is the default and preserves historical
  package behavior. It propagates exposure and event counts from the
  nearest non-empty interval within the same treatment group and emits a
  warning. This prevents undefined posterior draws, but it should be
  interpreted as a computational fallback rather than evidence about
  that later interval.
- `empty_interval = "prior"` leaves the interval with d\_{zj} = 0 and
  y\_{zj} = 0, so the posterior for that interval is exactly the
  specified Gamma prior.
- `empty_interval = "error"` stops the analysis when any treatment-arm
  interval has no exposed subjects.

The posterior density factorizes as

\pi(\boldsymbol{\lambda} \mid \mathcal{D}) = \prod_z \prod\_{j=1}^{J}
\pi(\lambda\_{zj} \mid d\_{zj}, y\_{zj}),

where each marginal factor is the Gamma distribution above. Posterior
predictive calculations integrate over this density rather than
conditioning on a single plug-in hazard estimate.

## 5. Predictive distribution for incomplete outcomes

At an interim analysis, subjects can be separated into three sets:

1.  enrolled subjects with complete endpoint information;
2.  enrolled subjects with partial follow-up; and
3.  future subjects not yet enrolled.

The first set contributes observed events and exposure to the posterior.
The second and third sets require prediction.

For an enrolled subject who is event-free through time u, a future event
time is drawn from the conditional piecewise-exponential distribution

\Pr(T \le t \mid T \> u) = \frac{F(t) - F(u)}{1 - F(u)}, \qquad t \> u,

where F(t) = 1 - S(t). Equivalently, if U \sim
\operatorname{Uniform}(0,1), then

T = F^{-1}\\F(u) + U\[1 - F(u)\]\\.

This is implemented by
[`pwe_impute()`](https://graemeleehickey.github.io/goldilocks/reference/pwe_impute.md).
For future subjects in the maximum-sample-size calculation,
[`pwe_sim()`](https://graemeleehickey.github.io/goldilocks/reference/pwe_sim.md)
draws unconditional event times from the same piecewise-exponential
model. In the package implementation, these future event-time
imputations do not require explicitly simulating future enrollment times
at the interim look; the number of future subjects is determined by
N\_{\max} - n\_\ell.

Let \mathcal{D}\_{\ell}^{\mathrm{obs}} denote the data observed at look
\ell, and let \mathcal{D}^{\mathrm{mis}} denote unobserved event times
and future censoring indicators. The posterior predictive density is

p(\mathcal{D}^{\mathrm{mis}} \mid \mathcal{D}\_{\ell}^{\mathrm{obs}}) =
\int p(\mathcal{D}^{\mathrm{mis}} \mid \boldsymbol{\lambda})
\pi(\boldsymbol{\lambda} \mid \mathcal{D}\_{\ell}^{\mathrm{obs}})
d\boldsymbol{\lambda}.

This density is the mathematical target, not an object that the package
evaluates in closed form. In the implementation, `posterior()` draws
\boldsymbol{\lambda}, `impute_data()` draws completed outcomes
conditional on those hazards, and `test_stop_success()` repeats the
impute-and-analyze cycle to approximate the predictive probability by
Monte Carlo simulation.

## 6. Interim decision algorithm

At interim look \ell, the package first estimates the posterior of the
hazard parameters from the currently observable data. It then uses Monte
Carlo integration to approximate P\_{n\_\ell} and P\_{\max,\ell}.

Let \psi(\mathcal{D}) be the final success indicator for a completed
analysis dataset:

\psi(\mathcal{D}) = I\\Q(\mathcal{D}) \> c\\,

where Q(\mathcal{D}) is the final analysis quantity and c is the success
threshold. In package notation, c is set by `prob_ha`. The final
analysis quantity Q(\mathcal{D}) depends on the analysis method:

| Design setting | `method` | Q(\mathcal{D}) | Supported alternatives |
|----|----|----|----|
| Two-arm randomized trial | `logrank` | 1-p(\mathcal{D}), where p(\mathcal{D}) is the traditional log-rank test P-value, with one-sided variants defined in Section 6.1 | `"less"`, `"greater"`, `"two.sided"` |
| Two-arm randomized trial | `cox` | 1-p(\mathcal{D}), where p(\mathcal{D}) is the traditional Wald-test P-value, with one-sided variants defined in Section 6.1 | `"less"`, `"greater"`, `"two.sided"` |
| Two-arm randomized trial | `riskdiff` | 1-p(\mathcal{D}), where p(\mathcal{D}) is the Wald-test P-value for the treatment-control event-risk difference | `"less"`, `"greater"`, `"two.sided"` |
| Two-arm randomized trial | `bayes-surv` | \Pr(\Delta \< h_0 \mid \mathcal{D}) or \Pr(\Delta \> h_0 \mid \mathcal{D}) | `"less"`, `"greater"` |
| Single-arm trial | `bayes-surv` | \Pr(p_1(\tau) \< h_0 \mid \mathcal{D}) or \Pr(p_1(\tau) \> h_0 \mid \mathcal{D}) | `"less"`, `"greater"` |
| Two-arm randomized trial | `bayes-bin` | \Pr(\Delta\_{\mathrm{bin}} \< h_0 \mid \mathcal{D}) or \Pr(\Delta\_{\mathrm{bin}} \> h_0 \mid \mathcal{D}) | `"less"`, `"greater"` |
| Single-arm trial | `bayes-bin` | \Pr(\pi_1 \< h_0 \mid \mathcal{D}) or \Pr(\pi_1 \> h_0 \mid \mathcal{D}) | `"less"`, `"greater"` |

For frequentist analyses, `prob_ha` is therefore a transformed P-value
threshold. For example, `prob_ha = 0.975` corresponds to a one-sided
\alpha = 0.025 rule. The value should be chosen during design
calibration to control the desired type I error rate across relevant
null scenarios.

### 5.1 Predictive probability at the current sample size

The current-sample-size predictive probability P\_{n\_\ell} is estimated
by:

1.  draw one set of hazard parameters from the interim posterior;
2.  impute remaining endpoint outcomes for enrolled subjects only;
3.  apply the prespecified final analysis rule to the completed data;
    and
4.  record whether that completed trial is successful.

Formally,

P\_{n\_\ell} = \operatorname{E}\\
\psi(\mathcal{D}\_{n\_\ell}^{\mathrm{comp}}) \mid
\mathcal{D}\_{\ell}^{\mathrm{obs}} \\,

where \mathcal{D}\_{n\_\ell}^{\mathrm{comp}} is the completed dataset
formed from the n\_\ell enrolled subjects after imputing their remaining
follow-up. Equivalently,

P\_{n\_\ell} = \int \psi(\mathcal{D}\_{n\_\ell}^{\mathrm{obs}},
\mathcal{D}\_{n\_\ell}^{\mathrm{mis}})
p(\mathcal{D}\_{n\_\ell}^{\mathrm{mis}} \mid
\mathcal{D}\_{\ell}^{\mathrm{obs}})
d\mathcal{D}\_{n\_\ell}^{\mathrm{mis}}.

Repeating the four-step procedure above for Monte Carlo replicate m =
1,\ldots,M, where M is set by `N_impute`, gives

\widehat{P}\_{n\_\ell} = \frac{1}{M}\sum\_{m=1}^{M} I\\\textrm{success
in replicate } m\\.

If

\widehat{P}\_{n\_\ell} \> S\_\ell,

accrual is stopped for expected success. Enrolled subjects are still
followed to the planned final analysis time.

### 5.2 Predictive probability at the maximum sample size

The maximum-sample-size predictive probability P\_{\max,\ell} is
estimated similarly, except that the completed trial includes both
currently enrolled subjects and future subjects required to reach
N\_{\max}. For each replicate, event times are imputed for the future
subjects, the completed dataset is analyzed, and success is recorded:

P\_{\max,\ell} = \operatorname{E}\\
\psi(\mathcal{D}\_{N\_{\max}}^{\mathrm{comp}}) \mid
\mathcal{D}\_{\ell}^{\mathrm{obs}} \\.

\widehat{P}\_{\max,\ell} = \frac{1}{M}\sum\_{m=1}^{M} I\\\textrm{success
at } N\_{\max} \textrm{ in replicate } m\\.

If

\widehat{P}\_{\max,\ell} \< F\_\ell,

the trial stops for futility. Otherwise, accrual continues to the next
interim look.

Thus the interim action at look \ell can be represented as

A\_\ell = \begin{cases} \textrm{stop accrual for expected success}, &
\widehat{P}\_{n\_\ell} \> S\_\ell,\\ \textrm{stop for futility}, &
\widehat{P}\_{\max,\ell} \< F\_\ell,\\ \textrm{continue accrual}, &
\textrm{otherwise}. \end{cases}

The first look satisfying either stopping condition defines the adaptive
stopping look,

L^\* = \inf\\\ell : \widehat{P}\_{n\_\ell} \> S\_\ell \textrm{ or }
\widehat{P}\_{\max,\ell} \< F\_\ell\\,

with L^\* = L + 1 if no interim stopping condition is met and the design
continues to N\_{\max}.

## 7. Final analysis

The final analysis is conducted after accrual has stopped and the
relevant follow-up has completed for the enrolled cohort, subject to the
handling of loss to follow-up described below. The final rule supplies
the binary success indicator used inside the predictive probability
calculations.

### 6.1 Frequentist final tests

For `method = "logrank"`, success is based on a log-rank test. For
`method = "cox"`, success is based on the Wald test from a Cox
proportional hazards regression. For `method = "riskdiff"`, success is
based on a Wald test for the treatment-control difference in binary
event risks at `end_of_study`. For these methods, `goldilocks` stores
1-p in `post_prob_ha`; this is not a posterior probability, but it puts
frequentist and Bayesian rules on a common “larger is stronger evidence”
scale. For example, a one-sided test at \alpha = 0.025 corresponds to
`prob_ha = 0.975`.

For the log-rank option, let Z\_{\mathrm{LR}} denote the signed log-rank
statistic, with positive values corresponding to excess events in the
control arm under the package convention. Let p\_{\mathrm{LR}} denote
the traditional two-sided log-rank test P-value. For the Cox option, let
\widehat{\eta} be the estimated log hazard ratio for treatment versus
control. The null value `h0` is on the log-hazard-ratio scale, so the
package uses

Z\_{\mathrm{Cox}} = \frac{\widehat{\eta} - h_0}
{\operatorname{se}(\widehat{\eta})}.

When `h0 = 0`, this is the usual hazard-ratio-equals-1 null. A
non-inferiority margin specified as a hazard ratio can be supplied as
`h0 = log(margin)`. Here, lower treatment hazard corresponds to
Z\_{\mathrm{Cox}} \< 0. Let p\_{\mathrm{Cox}} denote the two-sided
Wald-test P-value relative to `h0`. The package uses the following
method-specific definitions of Q(\mathcal{D}):

For the risk-difference option, let \widehat p_1 and \widehat p_0 be the
observed event proportions in the treatment and control arms, with
sample sizes n_1 and n_0. The estimated effect and its unpooled binomial
variance are

\widehat\Delta = \widehat p_1 - \widehat p_0,

U\_{\Delta} = \frac{\widehat p_1(1-\widehat p_1)}{n_1} + \frac{\widehat
p_0(1-\widehat p_0)}{n_0}.

The complete-data Wald statistic is

Z\_{\mathrm{RD}} = \frac{\widehat\Delta-h_0}{\sqrt{U\_{\Delta}}}.

| Method     | Alternative   | Q(\mathcal{D})                   |
|------------|---------------|----------------------------------|
| `logrank`  | `"less"`      | \Phi(Z\_{\mathrm{LR}})           |
| `logrank`  | `"greater"`   | 1 - \Phi(Z\_{\mathrm{LR}})       |
| `logrank`  | `"two.sided"` | 1 - p\_{\mathrm{LR}}             |
| `cox`      | `"less"`      | 1 - \Phi(Z\_{\mathrm{Cox}})      |
| `cox`      | `"greater"`   | \Phi(Z\_{\mathrm{Cox}})          |
| `cox`      | `"two.sided"` | 1 - p\_{\mathrm{Cox}}            |
| `riskdiff` | `"less"`      | 1 - \Phi(Z\_{\mathrm{RD}})       |
| `riskdiff` | `"greater"`   | \Phi(Z\_{\mathrm{RD}})           |
| `riskdiff` | `"two.sided"` | 1 - 2\Phi(-\|Z\_{\mathrm{RD}}\|) |

The one-sided directions differ between the log-rank rows and the
model-based rows because of the sign convention of the package’s
log-rank statistic.

The risk-difference analysis discards event-time information and
requires complete binary endpoint status. The returned `est_final` is
\widehat\Delta.

When a Cox or risk-difference final analysis uses multiple imputation,
the analysis is applied separately to each completed dataset. Let
\widehat{\theta}\_m and U_m be the scalar effect estimate and its
estimated variance from imputation m = 1,\ldots,M. For Cox regression
\widehat{\theta}\_m is the log hazard ratio; for risk difference it is
\widehat\Delta_m. Rubin’s scalar pooling rules give

\bar{\theta} = \frac{1}{M}\sum\_{m=1}^{M}\widehat{\theta}\_m, \qquad
\bar{U} = \frac{1}{M}\sum\_{m=1}^{M}U_m,

B = \frac{1}{M-1}\sum\_{m=1}^{M} (\widehat{\theta}\_m - \bar{\theta})^2,
\qquad T = \bar{U} + \left(1 + \frac{1}{M}\right)B.

The pooled Wald statistic is (\bar{\theta} - h_0) / \sqrt{T}. Its
P-value uses a t reference distribution with Rubin’s large-sample
degrees of freedom

\nu = (M-1)\left(1 + \frac{1}{r}\right)^2, \qquad r = \frac{(1 +
1/M)B}{\bar{U}}.

When B = 0, \nu = \infty and the reference distribution reduces to the
standard normal distribution. At least two imputations are therefore
required. The returned `est_final` is \bar{\theta}, while `post_prob_ha`
is 1-p from this pooled test with the direction determined by
`alternative`.

### 6.2 Bayesian survival final test

For `method = "bayes-surv"`, posterior hazard draws are mapped to
cumulative event probabilities at \tau. In a two-arm design the
treatment effect is

\Delta = p_1(\tau) - p_0(\tau),

where p_1(\tau) is the treatment-arm event probability and p_0(\tau) is
the control-arm event probability. The effect is on the event scale, not
the survival scale. For an adverse event, benefit usually means \Delta
\< 0.

Because p_a(\tau) = 1 - \exp\\-H_a(\tau)\\, posterior draws of \Delta
are obtained by transforming posterior draws of \boldsymbol{\lambda}\_1
and \boldsymbol{\lambda}\_0:

\Delta^{(b)} = \left\[1 - \exp\\-H_1^{(b)}(\tau)\\\right\] - \left\[1 -
\exp\\-H_0^{(b)}(\tau)\\\right\], \qquad b = 1,\ldots,B.

The Monte Carlo estimate of the posterior probability for
`alternative = "less"` is

\widehat{\Pr}(\Delta \< h_0 \mid \mathcal{D}) =
\frac{1}{B}\sum\_{b=1}^{B} I(\Delta^{(b)} \< h_0),

where B is set by `N_mcmc`.

With `alternative = "less"`, success is declared when

\Pr(\Delta \< h_0 \mid \mathcal{D}) \> \texttt{prob_ha}.

With `alternative = "greater"`, success is declared when

\Pr(\Delta \> h_0 \mid \mathcal{D}) \> \texttt{prob_ha}.

The Bayesian final test is one-sided in the package;
`alternative = "two.sided"` is not supported.

In a single-arm design there is no p_0(\tau). The estimand becomes
p_1(\tau), and h_0 is an external benchmark event probability. In
clinical-trial terminology this benchmark is often called a performance
goal (PG) or objective performance criterion (OPC). Consequently,
single-arm survival designs in `goldilocks` require
`method = "bayes-surv"`; complete binary single-arm designs can use
`method = "bayes-bin"`.

### 6.3 Bayesian binary final test

For `method = "bayes-bin"`, each analysis dataset is reduced to the
binary indicator of whether the endpoint has occurred by \tau. Subjects
with right-censored follow-up before \tau must be imputed or excluded
before this final test is applied, as described below.

Let x_z be the number of events and n_z the number of subjects in
treatment group z. With `bin_prior = c(a, b)`, the event probability in
arm z has posterior distribution

\pi_z \mid \mathcal{D} \sim \operatorname{Beta}(a + x_z, b + n_z - x_z).

In a two-arm design, the binary treatment effect is

\Delta\_{\mathrm{bin}} = \pi_1 - \pi_0,

the treatment-arm event probability minus the control-arm event
probability. For an adverse binary event, benefit usually means
\Delta\_{\mathrm{bin}} \< 0. With `alternative = "less"`, success is
declared when

\Pr(\Delta\_{\mathrm{bin}} \< h_0 \mid \mathcal{D}) \> \texttt{prob_ha}.

With `alternative = "greater"`, success is declared when

\Pr(\Delta\_{\mathrm{bin}} \> h_0 \mid \mathcal{D}) \> \texttt{prob_ha}.

In a single-arm design, the estimand is \pi_1 and `h0` is the external
benchmark event probability. Thus, `alternative = "less"` declares
success when

\Pr(\pi_1 \< h_0 \mid \mathcal{D}) \> \texttt{prob\\ha}.

The posterior probability can be computed in three ways. With
`bin_method = "mc"`, the package draws from the beta posterior directly.
With `bin_method = "normal"`, it uses a normal approximation to the
posterior mean or treatment-control difference. With
`bin_method = "quadrature"`, it uses numerical integration for the
two-arm posterior difference. The argument `N_mcmc` controls the number
of Monte Carlo beta draws only when `bin_method = "mc"`.

### 6.4 Loss to follow-up at the final analysis

Interim predictions impute outcomes that are not yet known. At the final
analysis, `imputed_final` controls whether subjects lost to follow-up
are also imputed.

If `imputed_final = TRUE`, Bayesian methods (`method = "bayes-surv"` or
`method = "bayes-bin"`) analyze each imputed completed dataset and
average the resulting posterior summaries. Cox regression and
risk-difference analyses instead pool completed-data scalar estimates
and variances using Rubin’s rules as described above; `N_impute` must be
at least two. Imputed final analyses remain unavailable for
`method = "logrank"` because no pooling rule is implemented for that
test.

If `imputed_final = FALSE`, the final analysis uses observed
right-censored data for methods that can handle censoring (`logrank`,
`cox`, and `bayes-surv`). For `riskdiff` and `bayes-bin`,
lost-to-follow-up subjects are excluded because these methods require
complete binary outcomes and have no mechanism for right-censored
observations. Rubin pooling applies to imputed Cox and risk-difference
final analyses; it does not alter the interim posterior-predictive
calculation, where each simulated completed trial is tested separately
before the success indicators are averaged.

The loss-to-follow-up mechanism in the simulator is non-informative.
Designs where dropout may depend on prognosis should be assessed with
sensitivity analyses outside the default data-generating mechanism.

## 8. Operating characteristics

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
- probability of reaching N\_{\max};
- mean and standard deviation of enrolled sample size; and
- probability of stopping for expected success but failing at the final
  analysis.

Let R = 1,\ldots,R\_{\max} index simulated trials under a scenario
\theta, where \theta denotes the data-generating parameters such as
control hazard, treatment hazard, accrual rate, loss-to-follow-up rate,
and follow-up duration. The trial-level random variables are:

| Symbol | Meaning                                                     |
|--------|-------------------------------------------------------------|
| N_R    | enrolled sample size in simulated trial R                   |
| E_R    | indicator that trial R stopped accrual for expected success |
| F_R    | indicator that trial R stopped for futility                 |
| Z_R    | final success indicator in trial R                          |

Let \Theta_0 denote the null parameter space, i.e. the set of
data-generating scenarios in which the treatment does not satisfy the
alternative hypothesis. For a two-arm superiority trial this includes
scenarios with no beneficial treatment effect; in practice, it should be
explored across plausible nuisance parameters such as control event
rates and accrual rates. The main operating characteristics are

\operatorname{Power}(\theta) = \Pr\_\theta(Z_R = 1, F_R = 0),

\operatorname{Type\\ I\\ error}(\theta_0) = \Pr\_{\theta_0}(Z_R = 1, F_R
= 0), \qquad \theta_0 \in \Theta_0,

\Pr\_\theta(\textrm{stop for expected success}) = \Pr\_\theta(E_R = 1),

\Pr\_\theta(\textrm{stop for futility}) = \Pr\_\theta(F_R = 1),

and

\operatorname{E}\_\theta(N_R), \qquad \operatorname{Var}\_\theta(N_R).

The `stop_and_fail` summary estimates

\Pr\_\theta(E_R = 1, Z_R = 0),

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
  method           = method,
  seed             = 12345)

summarise_sims(out$sims)
```

### 7.1 Visual diagnostics

The simulation plotting functions expose three different levels of the
design:

- [`plot_sim_ocs()`](https://graemeleehickey.github.io/goldilocks/reference/plot_sim_ocs.md)
  compares final success, stopping probabilities, and mean sample size
  across data-generating scenarios.
- [`plot_sim_stopping()`](https://graemeleehickey.github.io/goldilocks/reference/plot_sim_stopping.md)
  expands one scenario into marginal, conditional, or cumulative
  stopping summaries, or a count-based flowchart through successive
  looks.
- [`plot_sim_decisions()`](https://graemeleehickey.github.io/goldilocks/reference/plot_sim_decisions.md)
  examines the joint interim predictive probabilities and the decision
  thresholds at each look.

For operating-characteristic curves, first attach a numeric effect scale
to the scenario summary. The package does not infer this automatically
because the appropriate scale may be a hazard ratio, risk difference,
survival probability, or event probability depending on the analysis:

``` r

scenario_oc <- summarise_sims(list(
  "null" = null_sims$sims,
  "moderate" = moderate_sims$sims,
  "target" = target_sims$sims
))
scenario_oc$true_effect <- c(0, -0.10, -0.20)

plot_sim_ocs(
  scenario_oc,
  effect = "true_effect",
  xlab = "True treatment-control event-probability difference"
)
plot_sim_stopping(target_sims)
```

Decision maps require the optional simulation traces:

``` r

target_sims_traced <- update(target_sims, return_trace = TRUE)
plot_sim_stopping(target_sims_traced, type = "flowchart")
plot_sim_decisions(target_sims_traced)
```

Retaining traces does not change the trial summaries or random-number
path, but it increases the output size. Trace-recorded sample sizes also
let conditional, cumulative, and flowchart stopping views display
reached looks at which no trial stopped. Traces are therefore usually
most useful for selected scenarios after a broad
operating-characteristic grid has been screened.

Broglio et al. emphasize that type I error for this class of adaptive
design should be examined across the null space, not only at one
convenient null scenario. For time-to-event endpoints, the relevant null
space includes plausible control event rates and accrual rates. Accrual
rate is especially important because rapid enrollment can leave little
endpoint information available at interim looks, increasing the
uncertainty in both \widehat{P}\_{n\_\ell} and \widehat{P}\_{\max,\ell}.
The follow-up period after accrual stops also affects operating
characteristics because it determines how much additional information is
observed before the final analysis.

The Monte Carlo sizes `N_impute`, `N_mcmc`, and `N_trials` should be
chosen so that simulation error is small relative to the decision being
made. Small values are appropriate for examples and package tests, but
final design calibration usually requires many more trials and
imputations.

## 9. Threshold selection

The thresholds S\_\ell and F\_\ell may be constant across looks or may
vary by look. They interact with the final analysis threshold, number
and timing of looks, endpoint delay, accrual rate, loss to follow-up,
prior distribution, and maximum sample size.

If type I error is too high, possible remedies include increasing
`prob_ha`, increasing the expected-success thresholds, reducing the
number of looks, or altering follow-up requirements. If power is too
low, the maximum sample size, futility threshold, expected-success
threshold, or final analysis threshold may need reconsideration. Each
change should be rechecked under null and alternative scenarios.

The quantity `stop_and_fail` is particularly useful when tuning S\_\ell.
It estimates how often a trial stops accrual for expected success but
does not meet the final success criterion after follow-up is complete.
If this is too large, the expected-success threshold is usually too
permissive for the amount of uncertainty present at interim looks.

## 10. Relation to group-sequential designs

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

## 11. Package-specific scope

The time-to-event example in Broglio et al. (2014) uses a
Gamma-exponential prediction model. `goldilocks` extends this to a
piecewise-exponential model, so the hazard may change at prespecified
cut-points. This can be useful when there is a clinically plausible
early-risk period, but each additional interval adds parameters and can
make interim posteriors more diffuse.

The package also supports single-arm Bayesian Goldilocks designs by
replacing the concurrent control with an external benchmark h_0, often
referred to as a performance goal (PG) or objective performance
criterion (OPC). This is convenient for early-phase, rare-disease, or
proof-of-concept settings, but validity then depends on the benchmark
being transportable to the enrolled population.

For complete binary endpoints, `method = "bayes-bin"` replaces the
piecewise-exponential final analysis with a conjugate beta-binomial
analysis. The interim prediction machinery still imputes
not-yet-observed endpoint statuses from the event-time model, so users
should specify data-generating hazards that are clinically meaningful
for the binary endpoint time.

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
