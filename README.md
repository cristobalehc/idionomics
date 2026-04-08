# idionomics

[![R-hub](https://github.com/cristobalehc/idionomics/actions/workflows/rhub.yaml/badge.svg)](https://github.com/cristobalehc/idionomics/actions/workflows/rhub.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

**idionomics** is an R toolkit for **idionomic science** — a research philosophy that places the unit of the ensemble (individual/couple/group) at the center of analysis. Rather than assuming a common distribution, a similar enough process for each unit, and fitting a single model to the whole ensemble, idionomic methods model each unit separately, then aggregate upward if sensible. The group-level picture emerges from individual results, not the other way around, while explicitly evaluating whether aggregation is reasonable given the measured level of heterogeneity of effects. 

The package is built around intensive longitudinal data where each participant contributes a time series. It provides a pipeline from preprocessing through modeling to group-level summaries. 

---

## The idionomic science principle

Classical panel data methods (multilevel models, fixed-effects regression) estimate one set of parameters shared — or partially shared — across all units of an ensemble. If the focus of interest is the trajectory of individuals, this is only sensible under hard-to-meet assumptions, such as exchangeability and/or ergodicity. If these assumptions are not met, ensemble averages may systematically obscure individual differences: an average positive effect may coexist with a significant subset of individuals for whom the effect is negative, nonsignificant, or zero. 

Idionomic science inverts the order of operations:

1. **Unit first.** Fit a model to each person's time series independently, capturing that person's unique dynamics, residual structure, and effect sizes.
2. **Group later.** Aggregate the individual estimates with meta-analytic methods that explicitly represent and quantify heterogeneity across people (unsupervised clustering and other estimate-based methods are planned for future releases).

This preserves the individual's data structure, produces person-specific estimates that can be reported back to participants or explored as a basis for personalized intervention, and provides honest group-level summaries that distinguish "the average effect is X" from "most people show effect X", "the average effect is null but there are significant effects at both sides", and similar patterns that a single pooled estimate might drastically obscure.


---

## Installation

```r
# Install from a local source directory (development version)
install.packages("devtools")
devtools::install("/path/to/idionomics")

# Or install directly from GitHub:
devtools::install_github("cristobalehc/idionomics")
```

---

## Recommended analysis pipeline

```
i_screener()   →   pmstandardize()   →   i_detrender()   →   iarimax()/looping_machine()   →   i_pval() / sden_test()
```

1. **`i_screener()` [optional]** — pre-pipeline data quality filter. Removes or flags subjects with too few observations, insufficient raw variance, or repetitive responses before standardization. Should run on raw data, before `pmstandardize()`.
2. **`pmstandardize()` [optional]** — within-person z-scoring. 
3. **`i_detrender()` [optional]** — linear detrending. Removes the linear time trend within each subject and each variable independently, to reduce the differencing order auto.arima selects, among other uses.
4. **`iarimax()` / `looping_machine()`** — per-subject ARIMAX fitting and random-effects meta-analysis. `looping_machine()` extends this to three-variable directed loops (a→b, b→c, c→a).
5. **`i_pval()`** — attaches per-subject p-values based on ML-consistent degrees of freedom.
6. **`sden_test()`** — Sign Divergence / Equisyncratic Null test: a binomial test on the count of significant individual-level effects that are in the opposite direction of the pooled effect (Sign Divergence) or at both sides if pooled effect is not statistically significant (Equisyncratic Null test).


---

## Function reference

### `i_screener()` — Pre-pipeline data quality filter

```r
i_screener(df, cols, id_var,
         min_n_subject = 20,
         min_sd       = NULL,
         max_mode_pct = NULL,
         filter_type  = "joint",
         mode         = "filter",
         verbose      = FALSE)
```

Screens subjects for data quality on raw (unstandardised) data before it enters the pipeline. After `pmstandardize()`, all non-constant series have within-person variance = 1 by construction, making `iarimax()`'s `minvar` filter ineffective. Running `i_screener()` on raw data catches low-quality subjects at the right stage.

Three configurable criteria (all optional except `min_n_subject`):

| Criterion | What it catches | Default |
|---|---|---|
| `min_n_subject` | Subjects with too few observations | 20 |
| `min_sd` | Near-constant series (floor/ceiling, low range) | `NULL` (off) |
| `max_mode_pct` | "Stuck" responders (e.g. ≥ 80 % of responses identical) | `NULL` (off) |

`filter_type = "joint"` (default) excludes a subject if they fail any criterion on any variable — consistent with `iarimax()`'s AND filter. `filter_type = "per_column"` evaluates each variable independently.

`mode` controls the output format:
- `"filter"` — returns the dataframe with failing subjects removed (joint) or their failing column values set to `NA` (per_column).
- `"flag"` — appends a logical `pass_overall` column (joint) or `<col>_pass` columns (per_column).
- `"report"` — returns a per-subject quality summary table.

```r
library(idionomics)

set.seed(42)
panel <- do.call(rbind, lapply(1:9, function(id) {
  a <- rnorm(50)
  b <- 0.4 * a + rnorm(50)
  c <- 0.4 * b + rnorm(50)
  data.frame(id = as.character(id), time = seq_len(50),
             a = a, b = b, c = c,
             y = 0.5 * a + rnorm(50),
             stringsAsFactors = FALSE)
}))

# Subject 10: near-constant "a" — will be caught by i_screener(min_sd = 0.5)
s10 <- data.frame(id = "10", time = seq_len(50),
                  a = rep(3, 50), b = rnorm(50), c = rnorm(50),
                  y = rnorm(50), stringsAsFactors = FALSE)

# Subject 11: full positive loop (a -> b -> c -> a)
a11 <- rnorm(50)
b11 <- 0.6 * a11 + rnorm(50, sd = 0.5)
c11 <- 0.6 * b11 + rnorm(50, sd = 0.5)
s11 <- data.frame(id = "11", time = seq_len(50),
                  a = a11 + 0.4 * c11, b = b11, c = c11,
                  y = 0.5 * a11 + rnorm(50), stringsAsFactors = FALSE)

# Subject 12: negative a -> y effect
a12 <- rnorm(50)
s12 <- data.frame(id = "12", time = seq_len(50),
                  a = a12, b = 0.4 * a12 + rnorm(50),
                  c = rnorm(50), y = -0.5 * a12 + rnorm(50),
                  stringsAsFactors = FALSE)

panel <- rbind(panel, s10, s11, s12)

# Remove subjects with too few obs or low raw variance, before standardizing
panel_clean <- i_screener(panel, cols = c("a", "b", "c", "y"), id_var = "id",
                        min_n_subject = 20, min_sd = 0.5, verbose = TRUE)

# Inspect quality without committing to removal
report <- i_screener(panel, cols = c("a", "b", "c", "y"), id_var = "id",
                   min_n_subject = 20, min_sd = 0.5, max_mode_pct = 0.80,
                   mode = "report", verbose = TRUE)
print(report)

# Flag subjects for inspection, then decide
flagged <- i_screener(panel, cols = c("a", "b", "c", "y"), id_var = "id",
                    min_sd = 0.5, mode = "flag")
table(flagged$pass_overall)
```

---

### `pmstandardize()` — Within-person z-scoring

```r
pmstandardize(df, cols, id_var, verbose = FALSE, append = TRUE)
```

Computes `(x - person_mean) / person_sd` for each person × column combination. Output columns are named `<col>_psd`.

```r
# Standardize all four variables within each person
panel_std <- pmstandardize(panel_clean, cols = c("a", "b", "c", "y"), id_var = "id",
                          verbose = TRUE)
head(panel_std)
```

---

### `i_detrender()` — Within-person linear detrending

```r
i_detrender(df, cols, id_var, timevar,
            min_n_subject = 20, minvar = 0.01,
            verbose = FALSE, append = TRUE)
```

Fits `lm(col ~ time)` within each subject and appends the column with the residuals (`<col>_dt`). Subjects with too few observations, insufficient pre-detrend variance, or near-zero post-detrend variance receive `NA` — independently for each column.

```r
panel_dt <- i_detrender(panel_std, cols = c("a_psd", "b_psd", "c_psd", "y_psd"),
                        id_var = "id", timevar = "time", verbose = TRUE)
head(panel_dt)
```

---

### `iarimax()` — Core I-ARIMAX algorithm

```r
iarimax(dataframe, min_n_subject = 20, minvar = 0.01,
        y_series, x_series, focal_predictor = NULL,
        id_var, timevar, fixed_d = NULL,
        correlation_method = "pearson",
        keep_models = FALSE, verbose = FALSE)
```

Fits one `forecast::auto.arima()` model per subject, extracts coefficients via `broom::tidy()`, and pools the focal predictor's coefficients with `metafor::rma()`. The `fixed_d` argument optionally fixes the differencing order across all subjects to ensure coefficients are on the same scale (e.g., `fixed_d = 0` for levels, `fixed_d = 1` for changes); AR and MA orders are always selected automatically per subject.

```r
result <- iarimax(panel_dt,
                  y_series  = "y_psd_dt",
                  x_series  = "a_psd_dt",
                  id_var    = "id",
                  timevar   = "time",
                  verbose   = TRUE)

summary(result)   # prints subject counts, direction/significance counts, REMA
plot(result)      # caterpillar plot with RE-MA overlay
```

#### What the return value contains

| Field | Description |
|---|---|
| `$results_df` | Per-subject ARIMA orders, estimates, SEs, `n_valid`, `n_params`, `raw_cor` |
| `$meta_analysis` | `metafor::rma` object (or `NULL` if rma failed) |
| `$case_number_detail` | Subject counts: original, filtered, ARIMA-failed, analyzed |
| `$models` | Raw `Arima` objects (only if `keep_models = TRUE`) |

---

### `i_pval()` — Per-subject p-values

```r
i_pval(iarimax_object, feature = NULL)
```

Attaches a `pval_<feature>` column to `results_df` using the two-tailed t-distribution with ML-based degrees of freedom (`n_valid - n_params`).

```r
result_pval <- i_pval(result)
result_pval$results_df[, c("id", "estimate_a_psd_dt", "pval_a_psd_dt")]
```

---

### `sden_test()` — Sign Divergence / Equisyncratic Null test

```r
sden_test(iarimax_object, alpha_arimax = 0.05, alpha_binom = NULL,
          test = "auto", feature = NULL)
```

A binomial test on the count of individually significant effects. Two test variants:

- **ENT** (Equisyncratic Null Test): tests whether *any-direction* significant effects exceed chance. Selected automatically when the pooled REMA effect is non-significant.
- **SDT** (Sign Divergence Test): tests whether effects in the *counter-pooled* direction exceed chance. Selected automatically when the pooled REMA effect is significant. 

```r
sden <- sden_test(result)
summary(sden)

# Force ENT regardless of REMA
sden_ent <- sden_test(result, test = "ENT")
```

---

### `looping_machine()` — Directed loop detection

```r
looping_machine(dataframe, a_series, b_series, c_series, id_var, timevar,
                covariates = NULL, include_third_as_covariate = FALSE,
                min_n_subject = 20, minvar = 0.01, fixed_d = NULL,
                correlation_method = "pearson",
                alpha = 0.05, keep_models = FALSE, verbose = FALSE)
```

Fits three I-ARIMAX legs (a→b, b→c, c→a), applies `i_pval()` to each, and computes `Loop_positive_directed`: a 0/1 indicator that is 1 only when all three focal betas are positive *and* significant at `alpha`.

```r
loop_result <- looping_machine(panel_dt,
                               a_series = "a_psd_dt", b_series = "b_psd_dt",
                               c_series = "c_psd_dt",
                               id_var = "id", timevar = "time",
                               verbose = TRUE)

# Proportion of subjects with detected positive directed loop
mean(loop_result$loop_df$Loop_positive_directed, na.rm = TRUE)

# Per-leg I-ARIMAX results are also returned
summary(loop_result$iarimax_a_to_b)
```

---

## Full pipeline example

```r
library(idionomics)

set.seed(42)
panel <- do.call(rbind, lapply(1:9, function(id) {
  a <- rnorm(50)
  b <- 0.4 * a + rnorm(50)
  c <- 0.4 * b + rnorm(50)
  data.frame(
    id   = as.character(id),
    time = seq_len(50),
    a = a, b = b, c = c,
    y = 0.5 * a + rnorm(50),
    stringsAsFactors = FALSE
  )
}))

# Manually created subjects.
s10 <- data.frame(id = "10", time = seq_len(50),
                  a = rep(3, 50), b = rnorm(50), c = rnorm(50),
                  y = rnorm(50), stringsAsFactors = FALSE)
a11 <- rnorm(50)
b11 <- 0.6 * a11 + rnorm(50, sd = 0.5)
c11 <- 0.6 * b11 + rnorm(50, sd = 0.5)
s11 <- data.frame(id = "11", time = seq_len(50),
                  a = a11 + 0.4 * c11, b = b11, c = c11,
                  y = 0.5 * a11 + rnorm(50), stringsAsFactors = FALSE)
a12 <- rnorm(50)
s12 <- data.frame(id = "12", time = seq_len(50),
                  a = a12, b = 0.4 * a12 + rnorm(50),
                  c = rnorm(50), y = -0.5 * a12 + rnorm(50),
                  stringsAsFactors = FALSE)
panel <- rbind(panel, s10, s11, s12)

# Step 1: Quality screening on raw data (before standardization)
panel_clean <- i_screener(panel, cols = c("a", "b", "c", "y"), id_var = "id",
                        min_n_subject = 20, min_sd = 0.3, max_mode_pct = 0.80,
                        verbose = TRUE)

# Step 2: Within-person standardization
panel_std <- pmstandardize(panel_clean, cols = c("a", "b", "c", "y"), id_var = "id",
                          verbose = TRUE)

# Step 3: Linear detrending
panel_dt <- i_detrender(panel_std, cols = c("a_psd", "b_psd", "c_psd", "y_psd"),
                        id_var = "id", timevar = "time", verbose = TRUE)

# Step 4a: I-ARIMAX (single predictor)
result <- iarimax(panel_dt,
                  y_series = "y_psd_dt", x_series = "a_psd_dt",
                  id_var = "id", timevar = "time", verbose = TRUE)

summary(result)
plot(result)

# Step 4b: Directed loop detection
loop_result <- looping_machine(panel_dt,
                               a_series = "a_psd_dt", b_series = "b_psd_dt",
                               c_series = "c_psd_dt",
                               id_var = "id", timevar = "time",
                               verbose = TRUE)

mean(loop_result$loop_df$Loop_positive_directed, na.rm = TRUE)

# Step 5: Per-subject p-values
result_pval <- i_pval(result)

# Step 6: SDEN test
sden <- sden_test(result_pval)
summary(sden)
```

---

## Key dependencies

| Package | Purpose |
|---|---|
| `forecast` | `auto.arima()` — per-subject ARIMA model selection |
| `metafor` | `rma()` — random-effects meta-analysis |
| `broom` | `tidy()` — coefficient extraction from Arima objects |
| `ggplot2`, `forcats` | Caterpillar plot |
| `dplyr`, `tidyr`, `tibble`, `rlang` | Data manipulation |

---

## Disclaimer

**idionomics** is free, open-source software provided "as-is", without warranty of any kind — express or implied — under the [MIT License](LICENSE.md). The authors are not liable for any damages or losses arising from the use of this software.

Statistical software can produce results that are technically valid but analytically inappropriate for a given context. Users are encouraged to review the methods and code, inspect their data, and exercise independent statistical judgment before reporting findings.
