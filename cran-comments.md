## Submission notes

The win-builder R-oldrelease check reported one incoming-feasibility NOTE.
The reported issues have been addressed:

* `DESCRIPTION` now uses `Broglio and colleagues (2014)`, avoiding the words
  `et al.` flagged by the spelling check and retaining the citation wording
  used in the previous release.
* Both links to the Wilber trial report in the ThermoCool AF vignette now use
  its DOI, `10.1001/jama.2009.2029`, instead of the PubMed URL that returned
  HTTP 203.

The revised source archive passes `R CMD check --as-cran` on macOS with
R 4.6.1: 0 errors, 0 warnings, and 0 notes.

This is goldilocks 1.0.0, the first stable major release and an update to
CRAN version 0.6.0. The main changes are:

* observed-interim evaluation through `evaluate_interim()`;
* restricted mean survival time and Farrington-Manning analyses;
* an optional immediate-success stopping boundary;
* arm-specific priors, separate generation and analysis partitions, and
  corrected predictive/final analysis-prior handling;
* independent exponential dropout and improved simulation robustness;
* Monte Carlo uncertainty and calendar-time operating-characteristic summaries;
* additional statistical validation and expanded worked examples; and
* a Zenodo badge and smaller HTML vignettes with shared local diagram resources.

NEWS.md starts with migration guidance from 0.6.0. It identifies changes to
dropout, empty-interval handling, Bayesian survival predictions, Monte Carlo
defaults, method names, and direct `sim_comp_data()` calls. Users should
recalibrate existing designs after updating. The maintainer and license are
unchanged.

## Test environments

* Local macOS Tahoe 26.6.2, Apple silicon, R 4.6.1 (current release), using
  `R CMD build --compact-vignettes=gs+qpdf` and `R CMD check --as-cran`.
  R's standard compiler settings were used with `R_MAKEVARS_USER=/dev/null`.
  HTML Tidy 5.8.0 and V8 were available for complete HTML/manual validation.
* The implementation at commit 651944884215ee7edd9d28cea854856eeffa73b7 passed
  GitHub Actions on Windows release, macOS release, Ubuntu release,
  Ubuntu oldrel-1, and Ubuntu devel:
  https://github.com/graemeleehickey/goldilocks/actions/runs/34261133547
  These CI results precede the release metadata and vignette packaging changes;
  the final source archive was built and checked locally with R 4.6.1.

R-devel is not installed locally; the final archive was checked with the
current R release, supplemented by the implementation's R-devel CI result.

## R CMD check results

0 errors | 0 warnings | 0 notes

* All 1,506 CRAN test expectations pass; eight longer tests are skipped on CRAN.
* The full installed-package test suite passes all 1,527 expectations with no
  failures, warnings, or skips, including statistical calibration and parallel
  worker checks.
* All examples and 14 vignettes build and rebuild successfully. PDF and HTML
  manual checks, including HTML validation and math rendering, pass.
* The source archive is about 2.5 MB. Installed documentation is about 4.1 MB;
  the diagram JavaScript is shared rather than embedded repeatedly. All local
  vignette resource links were checked, and diagram rendering was verified in
  a browser.

## URL checks

The CRAN incoming URL check passes. A separate `urlchecker` run encountered
publisher access errors for these registered DOIs:

* 10.1161/CIRCHEARTFAILURE.119.005879 (HTTP 403; bibliographic record verified
  at https://pubmed.ncbi.nlm.nih.gov/31722536/).
* 10.1200/JCO.2014.55.2208 (HTTP 403; article verified
  at https://pmc.ncbi.nlm.nih.gov/articles/PMC4105489/).
* 10.1201/EBK1439825488 (HTTP 503; Crossref bibliographic record verified).

Zenodo and PubMed also returned intermittent timeouts/nonstandard responses.
The requested Zenodo DOI and SVG badge both resolved with HTTP 200 during
initial direct verification; later requests intermittently timed out.
The DOI identifies the archived 0.6.0
release and is retained as supplied; it does not claim to archive 1.0.0.

## Reverse dependencies

No reverse dependencies (Depends, Imports, LinkingTo, Suggests, or Enhances)
were found in the CRAN package database on 2026-09-08. CRAN's published checks
for version 0.6.0 were also all OK on that date.
