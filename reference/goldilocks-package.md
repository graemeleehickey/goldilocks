# goldilocks: Adaptive Trial Designs for Survival and Binary Endpoints

Implements Goldilocks adaptive trial designs for time-to-event and
fixed-time binary endpoints. Outcomes are generated with a piecewise
exponential model, with conjugate Gamma priors used for predictive
imputation. Final analyses may use log-rank or Cox tests, Bayesian
piecewise-exponential inference, frequentist risk differences, or
Bayesian beta-binomial inference. The method closely follows Broglio and
colleagues
[doi:10.1080/10543406.2014.888569](https://doi.org/10.1080/10543406.2014.888569)
and supports simulation of design operating characteristics.

## See also

Useful links:

- <https://graemeleehickey.github.io/goldilocks/>

- <https://github.com/graemeleehickey/goldilocks>

- <https://CRAN.R-project.org/package=goldilocks>

- Report bugs at <https://github.com/graemeleehickey/goldilocks/issues>

## Author

**Maintainer**: Graeme L. Hickey <graemeleehickey@gmail.com>
([ORCID](https://orcid.org/0000-0002-4989-0054))

Authors:

- Ying Wan <ying.wan@bd.com>

- Thevaa Chandereng <tc3123@cumc.columbia.edu>
  ([ORCID](https://orcid.org/0000-0003-4078-9176)) (bayesDP code as a
  template)

Other contributors:

- Becton, Dickinson and Company \[copyright holder\]

- Tim Kacprowski <t.kacprowski@tu-braunschweig.de> (For code from
  fastlogrank R package.) \[contributor\]
