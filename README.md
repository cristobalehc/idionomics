# idionomicsv2

**idionomicsv2** is an R toolkit for **idionomic science** — a research philosophy that places the individual at the center of analysis. Rather than fitting a single model to the group and treating between-person averages as the primary finding, idionomic methods model each person separately, then aggregate upward. The group-level picture emerges from individual results, not the other way around.

The package is built around intensive longitudinal data (ILD): experience-sampling, ecological momentary assessment, daily diary, wearable sensor streams, and similar designs where each participant contributes a time series. It provides the full pipeline from preprocessing through modeling to population-level inference.

---

## The idionomic science principle

Classical panel data methods (multilevel models, fixed-effects regression) estimate one set of parameters shared — or partially shared — across all subjects. If the focus of interest is the trajectory of individuals, this is only sensible under hard to meet assumptions such as exchangeability and/or ergodicity. If these assumptions are not met, ensemble averages may systematically obscure individual differences: an average positive effect may coexist with a significant subset of individuals for whom the effect is negative, nonsignificant, or zero. 

Idionomic science inverts the order of operations:

1. **Individual first.** Fit a model to each person's time series independently, capturing that person's unique dynamics, autocorrelation structure, and effect sizes.
2. **Group later.** Aggregate the individual estimates with meta-analytic, unsupervised clustering or estimate-based methods that explicitly represent and quantify heterogeneity across people.

This preserves the individual's data structure, produces person-specific estimates that can be reported back to participants or used for personalized intervention, and provides honest group-level summaries that distinguish "the average effect is X" from "most people show effect X", "the average effect is null but there are significant effects at both sides", and similar patterns that a single pooled estimate would erase.


---

## Installation

```r
# Install from a local source directory (development version)
install.packages("devtools")
devtools::install("/path/to/idionomicsv2")

# Or install directly from GitHub:
devtools::install_github("cristobalehc/idionomics")
```

---

## Recommended analysis pipeline

```
pmstandardize()   →   i_detrender()   →   iarimax()   →   i_pval() / sden_test()
```

1. **`pmstandardize()`** — within-person z-scoring. Removes between-person differences in mean and variance so coefficients are comparable across subjects.
2. **`i_detrender()`** — linear detrending. Removes the linear time trend within each subject and each variable independently, so that a spurious trend does not inflate the association estimate.
3. **`iarimax()`** — per-subject ARIMAX fitting and random-effects meta-analysis.
4. **`i_pval()`** — attaches per-subject p-values based on ML-consistent degrees of freedom.
5. **`sden_test()`** — Sign Divergence / Equisyncratic Null test: a binomial test on the count of significant individual-level effects.

> **Order matters.** Always standardize before detrending. If you detrend first and then standardize, near-zero residuals divided by a near-zero SD get rescaled to apparent unit variance, bypassing `iarimax()`'s variance filter and inflating coefficients.

---

## Function reference

### `pmstandardize()` — Within-person z-scoring

```r
pmstandardize(df, cols, idvar, verbose = FALSE, append = TRUE)
```

Computes `(x - person_mean) / person_sd` for each person × column combination. Output columns are named `<col>_psd`.

```r
library(idionomicsv2)

set.seed(1)
panel <- do.call(rbind, lapply(1:4, function(id) {
  data.frame(id = as.character(id), time = seq_len(25),
             x = rnorm(25, mean = id), y = rnorm(25),
             stringsAsFactors = FALSE)
}))

# Standardize x and y within each person
panel_std <- pmstandardize(panel, cols = c("x", "y"), idvar = "id")
head(panel_std)
```

---

### `i_detrender()` — Within-person linear detrending

```r
i_detrender(df, cols, idvar, timevar,
            min_n_subject = 20, minvar = 0.01,
            verbose = FALSE, append = TRUE)
```

Fits `lm(col ~ time)` within each subject and replaces the column with the residuals (`<col>_DT`). Subjects with too few observations, insufficient pre-detrend variance, or near-zero post-detrend variance receive `NA` — independently for each column.

```r
panel_dt <- i_detrender(panel_std, cols = c("x_psd", "y_psd"),
                        idvar = "id", timevar = "time")
head(panel_dt)
```

---

### `iarimax()` — Core I-ARIMAX algorithm

```r
iarimax(dataframe, y_series, x_series, id_var, timevar,
        focal_predictor = NULL,
        min_n_subject = 20, minvar = 0.01,
        correlation_method = "pearson",
        keep_models = FALSE, verbose = FALSE)
```

Fits one `auto.arima()` model per subject, extracts coefficients via `broom::tidy()`, and pools the focal predictor's coefficients with `metafor::rma()`.

```r
result <- iarimax(panel_dt,
                  y_series  = "y_psd_DT",
                  x_series  = "x_psd_DT",
                  id_var    = "id",
                  timevar   = "time")

summary(result)   # prints subject counts, direction/significance counts, REMA
plot(result)      # caterpillar plot with RE-MA overlay
```

#### What the return value contains

| Field | Description |
|---|---|
| `$results_df` | Per-subject ARIMA orders, estimates, SEs, `n_valid`, `n_params`, `raw_cor` |
| `$meta_analysis` | `metafor::rma` object (or `NULL` if fewer than 2 valid models) |
| `$case_number_detail` | Subject counts: original, filtered, ARIMA-failed, analyzed |
| `$models` | Raw `Arima` objects (only if `keep_models = TRUE`) |

---

### `i_pval()` — Per-subject p-values

```r
i_pval(iarimax_object, feature = NULL)
```

Attaches a `pval_<feature>` column to `results_df` using the two-tailed t-distribution with ML-based degrees of freedom (`n_valid - n_params`). P-values differ from `lm()` because ARIMA uses `sigma² = SSR/n`, not `SSR/(n-k)`.

```r
result_pval <- i_pval(result)
result_pval$results_df[, c("id", "estimate_x_psd_DT", "pval_x_psd_DT")]
```

---

### `sden_test()` — Sign Divergence / Equisyncratic Null test

```r
sden_test(iarimax_object, alpha_arimax = 0.05, alpha_binom = NULL,
          test = "auto", feature = NULL)
```

A binomial test on the count of individually significant effects. Two test variants:

- **ENT** (Equisyncratic Null Test): tests whether *any-direction* significant effects exceed chance. Selected automatically when the pooled REMA effect is non-significant.
- **SDT** (Sign Divergence Test): tests whether effects in the *counter-pooled* direction exceed chance. Selected automatically when the pooled REMA effect is significant, as a sensitivity check for heterogeneity.

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
                min_n_subject = 20, minvar = 0.01,
                correlation_method = "pearson",
                alpha = 0.05, keep_models = FALSE, verbose = FALSE)
```

Fits three I-ARIMAX legs (a→b, b→c, c→a), applies `i_pval()` to each, and computes `Loop_positive_directed`: a 0/1 indicator that is 1 only when all three focal betas are positive *and* significant at `alpha`.

```r
set.seed(7)
panel3 <- do.call(rbind, lapply(1:6, function(id) {
  a <- rnorm(30)
  b <- 0.4 * a + rnorm(30)
  c <- 0.4 * b + rnorm(30)
  data.frame(id = as.character(id), time = seq_len(30),
             a = a, b = b, c = c, stringsAsFactors = FALSE)
}))

loop_result <- looping_machine(panel3,
                               a_series = "a", b_series = "b", c_series = "c",
                               id_var = "id", timevar = "time")

# Proportion of subjects with detected positive directed loop
mean(loop_result$loop_df$Loop_positive_directed, na.rm = TRUE)

# Per-leg I-ARIMAX results are also returned
summary(loop_result$iarimax_a_to_b)
```

---

## Full pipeline example

```r
library(idionomicsv2)

set.seed(42)
panel <- do.call(rbind, lapply(1:10, function(id) {
  x <- rnorm(50)
  data.frame(
    id   = as.character(id),
    time = seq_len(50),
    x    = x,
    y    = 0.4 * x + rnorm(50),
    stringsAsFactors = FALSE
  )
}))

# Step 1: Within-person standardization
panel_std <- pmstandardize(panel, cols = c("x", "y"), idvar = "id")

# Step 2: Linear detrending
panel_dt <- i_detrender(panel_std, cols = c("x_psd", "y_psd"),
                        idvar = "id", timevar = "time")

# Step 3: I-ARIMAX
result <- iarimax(panel_dt,
                  y_series = "y_psd_DT", x_series = "x_psd_DT",
                  id_var = "id", timevar = "time")

# Step 4: Summary and plot
summary(result)
plot(result, y_series_name = "Mood", x_series_name = "Stress")

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

## Package status

Current version: `0.0.0.9000` (development)

`R CMD check` result: `0 errors | 0 warnings | 0 notes`
