# goldilocks

The goal of `goldilocks` is to implement the Goldilocks Bayesian
adaptive design proposed by Broglio et al. (2014), for one- and two-arm
trials. Outcomes are generated with an underlying piecewise-exponential
event-time model. Final analyses may retain the time-to-event outcome
(`method = "logrank"`, `"cox"`, or `"bayes-surv"`) or reduce complete
follow-up to event status at a fixed endpoint time
(`method = "riskdiff"` or `"bayes-bin"`).

The method can be used for a confirmatory trial to select a sample size
based on accumulating data. During accrual, predictive probabilities are
used to determine whether the current sample size is sufficient, whether
continuing accrual would be futile, or whether enrollment should
continue. The algorithm explicitly accounts for completion of planned
follow-up before the primary analysis. Broglio et al. (2014) refer to
this as a Goldilocks trial design, as it is constantly asking the
question, **“Is the sample size too big, too small, or just right?”**

## References

Broglio KR, Connor JT, Berry SM. Not too big, not too small: a
Goldilocks approach to sample size selection. *Journal of
Biopharmaceutical Statistics*, 2014; **24(3)**: 685–705.
