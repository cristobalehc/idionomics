# Tests for i_pval()
# Layer 1: formula verification against manual t-test calculation (fast, no auto.arima)
# Layer 2: integration comparison against lm() for white-noise data (skip_on_cran)

# ── Helper: minimal fake iarimax_results ──────────────────────────────────────
# Builds a fake object without running iarimax(), so Layer 1 tests are fast
# and not sensitive to auto.arima behaviour.

make_fake_iarimax <- function(include_na = FALSE) {
  results_df <- data.frame(
    id             = c("1", "2", "3"),
    estimate_x     = c(0.5, -0.3, 0.8),
    "std.error_x"  = c(0.1,  0.15, 0.2),
    n_valid        = c(25L,  25L,  25L),
    n_params       = c(2L,   2L,   2L),
    stringsAsFactors = FALSE,
    check.names    = FALSE
  )

  if (include_na) {
    # Simulate a subject whose ARIMA fit failed (all estimate columns are NA)
    results_df[2, "estimate_x"]    <- NA
    results_df[2, "std.error_x"]   <- NA
  }

  obj <- list(
    results_df           = results_df,
    meta_analysis        = NULL,
    error_arimax_skipped = character(0),
    models               = NULL
  )
  class(obj) <- c("iarimax_results", "list")
  attr(obj, "focal_predictor") <- "x"
  attr(obj, "id_var")          <- "id"
  attr(obj, "timevar")         <- "time"
  obj
}

# ── Helper: white-noise panel for lm comparison ───────────────────────────────
# Pure white noise residuals so auto.arima reliably selects ARIMA(0,0,0).

make_wn_panel <- function(n_subjects = 4, n_obs = 60, seed = 7) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(n_subjects), function(id) {
    x <- rnorm(n_obs)
    y <- 0.5 * x + rnorm(n_obs)
    data.frame(id = as.character(id), time = seq_len(n_obs),
               x = x, y = y, stringsAsFactors = FALSE)
  }))
}


# ══════════════════════════════════════════════════════════════════════════════
# Layer 1 — Formula verification (fast, no CRAN skip needed)
# ══════════════════════════════════════════════════════════════════════════════

test_that("i_pval p-values match manual two-tailed t-test formula", {
  fake     <- make_fake_iarimax()
  result   <- i_pval(fake)
  df       <- fake$results_df

  expected <- 2 * stats::pt(
    -abs(df$estimate_x / df[["std.error_x"]]),
    df = df$n_valid - df$n_params
  )

  expect_equal(result$results_df$pval_x, expected, tolerance = 1e-12)
})

test_that("p-values are in (0, 1] for valid estimates", {
  result <- i_pval(make_fake_iarimax())
  pvals  <- result$results_df$pval_x
  expect_true(all(pvals > 0 & pvals <= 1, na.rm = TRUE))
})

test_that("positive and negative estimates of equal magnitude give same p-value", {
  fake <- make_fake_iarimax()
  # subjects 1 and 3 have estimates 0.5 and 0.8 — make them symmetric
  fake$results_df$estimate_x    <- c(0.5, -0.5, 0.5)
  fake$results_df[["std.error_x"]] <- c(0.1,  0.1,  0.1)

  result <- i_pval(fake)
  expect_equal(result$results_df$pval_x[1],
               result$results_df$pval_x[2],
               tolerance = 1e-12)
})

test_that("NA estimate produces NA p-value, valid subjects unaffected", {
  fake   <- make_fake_iarimax(include_na = TRUE)
  result <- i_pval(fake)
  pvals  <- result$results_df$pval_x

  expect_true(is.na(pvals[2]))             # failed subject → NA
  expect_false(is.na(pvals[1]))            # valid subjects unaffected
  expect_false(is.na(pvals[3]))
})

test_that("i_pval only adds the pval column and does not modify other fields", {
  fake   <- make_fake_iarimax()
  result <- i_pval(fake)

  # pval column added
  expect_true("pval_x" %in% names(result$results_df))

  # estimates and std.errors unchanged
  expect_equal(result$results_df$estimate_x,      fake$results_df$estimate_x)
  expect_equal(result$results_df[["std.error_x"]], fake$results_df[["std.error_x"]])
  expect_equal(result$results_df$n_valid,         fake$results_df$n_valid)
  expect_equal(result$results_df$n_params,        fake$results_df$n_params)

  # other object fields untouched
  expect_null(result$meta_analysis)
  expect_null(result$models)
  expect_equal(result$error_arimax_skipped, character(0))
})

test_that("i_pval works for a non-focal predictor in a multi-predictor object", {
  # Build a fake object with two predictors
  fake <- make_fake_iarimax()
  fake$results_df$estimate_z     <- c(1.0, -0.5, 0.3)
  fake$results_df[["std.error_z"]] <- c(0.2,  0.1,  0.4)

  result <- i_pval(fake, feature = "z")

  expected <- 2 * stats::pt(
    -abs(fake$results_df$estimate_z / fake$results_df[["std.error_z"]]),
    df = fake$results_df$n_valid - fake$results_df$n_params
  )
  expect_equal(result$results_df$pval_z, expected, tolerance = 1e-12)
  # focal predictor pval column not added
  expect_false("pval_x" %in% names(result$results_df))
})

test_that("i_pval errors on wrong object class", {
  expect_error(i_pval(list(results_df = data.frame())),
               "iarimax_results")
})

test_that("i_pval errors when feature column is missing", {
  fake <- make_fake_iarimax()
  expect_error(i_pval(fake, feature = "nonexistent"),
               "not found")
})

test_that("i_pval returns NA (not NaN) and warns when df <= 0", {
  # n_valid == n_params => df = 0 for subject 1
  fake <- make_fake_iarimax()
  fake$results_df$n_valid  <- c(2L, 25L, 25L)
  fake$results_df$n_params <- c(2L,  2L,  2L)

  expect_warning(result <- i_pval(fake), regexp = "degrees of freedom")

  pvals <- result$results_df$pval_x
  expect_true(is.na(pvals[1]))    # df = 0  -> NA, not NaN
  expect_false(is.na(pvals[2]))   # df = 23 -> valid p-value
  expect_false(is.na(pvals[3]))
})

test_that("each subject's p-value corresponds to its own estimate and SE, not a neighbour's", {
  # Three subjects with deliberately distinct |t| = 20, 3, 0.5 so p-values are
  # well-separated. We verify each position against its own formula individually.
  fake <- make_fake_iarimax()
  fake$results_df$estimate_x    <- c( 2.0,  0.3, 0.05)
  fake$results_df[["std.error_x"]] <- c( 0.1,  0.1,  0.1)
  fake$results_df$n_valid       <- c(25L, 25L, 25L)
  fake$results_df$n_params      <- c( 2L,  2L,  2L)

  result <- i_pval(fake)
  pvals  <- result$results_df$pval_x

  expect_equal(pvals[1], 2 * stats::pt(-abs(2.0  / 0.1), df = 23), tolerance = 1e-12)
  expect_equal(pvals[2], 2 * stats::pt(-abs(0.3  / 0.1), df = 23), tolerance = 1e-12)
  expect_equal(pvals[3], 2 * stats::pt(-abs(0.05 / 0.1), df = 23), tolerance = 1e-12)
})


# ══════════════════════════════════════════════════════════════════════════════
# Layer 2 — LM comparison (integration, skip on CRAN)
# ══════════════════════════════════════════════════════════════════════════════

skip_on_cran()

# Note: lm() is NOT used as ground truth here. ARIMA uses ML estimation for
# sigma^2 (SSR/n) while lm() uses OLS (SSR/(n-k)), producing systematically
# different standard errors and p-values even for ARIMA(0,0,0). The correct
# ground truth is the formula itself, applied to real iarimax output.

panel  <- make_wn_panel()
result <- iarimax(dataframe = panel, y_series = "y", x_series = "x",
                  id_var = "id", timevar = "time")

test_that("i_pval formula holds on real iarimax output", {
  result_pval <- i_pval(result)
  df          <- result_pval$results_df

  valid    <- !is.na(df$estimate_x) & !is.na(df[["std.error_x"]])
  expected <- 2 * stats::pt(
    -abs(df$estimate_x[valid] / df[["std.error_x"]][valid]),
    df = df$n_valid[valid] - df$n_params[valid]
  )

  expect_equal(df$pval_x[valid], expected, tolerance = 1e-12)
})

test_that("i_pval p-values are in (0, 1] on real iarimax output", {
  result_pval <- i_pval(result)
  pvals       <- result_pval$results_df$pval_x
  expect_true(all(pvals > 0 & pvals <= 1, na.rm = TRUE))
})

test_that("larger |t-statistic| produces smaller p-value on real iarimax output", {
  result_pval <- i_pval(result)
  df          <- result_pval$results_df

  valid  <- !is.na(df$estimate_x) & !is.na(df[["std.error_x"]])
  t_abs  <- abs(df$estimate_x[valid] / df[["std.error_x"]][valid])
  pvals  <- df$pval_x[valid]

  # rank correlation between |t| and p-value should be strongly negative
  expect_true(cor(t_abs, pvals, method = "spearman") < -0.9)
})
