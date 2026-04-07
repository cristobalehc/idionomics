# idionomics 0.1.0

Initial CRAN release.

## Core algorithm

* `iarimax()` — per-subject ARIMAX fitting via `forecast::auto.arima()` followed
  by random-effects meta-analysis via `metafor::rma()`. Supports multi-predictor
  models, optional `fixed_d` for cross-subject comparability, and `keep_models`
  for retaining raw model objects.

## Preprocessing

* `i_screener()` — pre-pipeline data quality screening on raw data. Three
  criteria: minimum observations (`min_n_subject`), minimum within-person SD
  (`min_sd`), and maximum modal response proportion (`max_mode_pct`). Supports
  `"filter"`, `"flag"`, and `"report"` output modes.

* `pmstandardize()` — within-person z-scoring (person-mean centering and
  person-SD scaling).

* `i_detrender()` — within-person linear detrending via `lm(col ~ timevar)`.
  Per-column filtering with pre- and post-detrend variance guards.

## Inference

* `i_pval()` — per-subject p-values using the two-tailed t-distribution with
  ML-based degrees of freedom (`n_valid - n_params`).

* `sden_test()` — Sign Divergence Test (SDT) and Equisyncratic Null Test (ENT)
  with automatic selection based on the pooled REMA p-value.

## Loop detection

* `looping_machine()` — directed loop detection across three variables. Fits
  three `iarimax()` legs, applies `i_pval()`, and computes
  `Loop_positive_directed`.

## S3 methods

* `summary.iarimax_results()` — subject counts, direction/significance counts,
  REMA estimates, and heterogeneity statistics.

* `plot.iarimax_results()` — caterpillar plot with per-subject confidence
  intervals and REMA band overlay.

* `summary.sden_results()` — test type, hypothesis, and binomial test results.
