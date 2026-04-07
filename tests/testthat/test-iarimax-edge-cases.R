# Edge cases and boundary conditions.
# Each test calls iarimax() directly and has its own skip_on_cran().

# ── Filtering: min_n_subject ──────────────────────────────────────────────────

test_that("subject below min_n_subject is absent from results_df", {
  skip_on_cran()
  base  <- make_panel(n_subjects = 2, n_obs = 25)
  short <- data.frame(id = "short", time = seq_len(10),
                      x = rnorm(10), y = rnorm(10),
                      stringsAsFactors = FALSE)
  panel <- rbind(base, short)

  res <- iarimax(dataframe = panel, y_series = "y", x_series = "x",
                 id_var = "id", timevar = "time", min_n_subject = 20)

  expect_false("short" %in% res$results_df$id)
  expect_true(all(c("1", "2") %in% res$results_df$id))
})

test_that("subject with exactly min_n_subject observations is included", {
  skip_on_cran()
  base    <- make_panel(n_subjects = 2, n_obs = 25)
  on_edge <- data.frame(id = "edge", time = seq_len(20),
                        x = rnorm(20), y = rnorm(20),
                        stringsAsFactors = FALSE)
  panel <- rbind(base, on_edge)

  res <- iarimax(dataframe = panel, y_series = "y", x_series = "x",
                 id_var = "id", timevar = "time", min_n_subject = 20)

  expect_true("edge" %in% res$results_df$id)
})

# ── Filtering: minvar ─────────────────────────────────────────────────────────

test_that("subject with constant y is absent from results_df", {
  skip_on_cran()
  base  <- make_panel(n_subjects = 2, n_obs = 25)
  flat  <- data.frame(id = "flat_y", time = seq_len(25),
                      x = rnorm(25), y = rep(5, 25),
                      stringsAsFactors = FALSE)
  panel <- rbind(base, flat)

  res <- iarimax(dataframe = panel, y_series = "y", x_series = "x",
                 id_var = "id", timevar = "time")

  expect_false("flat_y" %in% res$results_df$id)
})

test_that("subject with constant x is absent from results_df", {
  skip_on_cran()
  base  <- make_panel(n_subjects = 2, n_obs = 25)
  flat  <- data.frame(id = "flat_x", time = seq_len(25),
                      x = rep(3, 25), y = rnorm(25),
                      stringsAsFactors = FALSE)
  panel <- rbind(base, flat)

  res <- iarimax(dataframe = panel, y_series = "y", x_series = "x",
                 id_var = "id", timevar = "time")

  expect_false("flat_x" %in% res$results_df$id)
})

# ── Filtering: NAs reduce effective n ────────────────────────────────────────

test_that("subject with many NA y rows can fall below min_n_subject", {
  skip_on_cran()
  # 25 rows but 12 have NA y -> only 13 complete cases -> filtered at min_n=20
  base  <- make_panel(n_subjects = 2, n_obs = 25)
  set.seed(1)
  sparse <- data.frame(id = "sparse", time = seq_len(25),
                       x = rnorm(25), y = rnorm(25),
                       stringsAsFactors = FALSE)
  sparse$y[1:12] <- NA
  panel <- rbind(base, sparse)

  res <- iarimax(dataframe = panel, y_series = "y", x_series = "x",
                 id_var = "id", timevar = "time", min_n_subject = 20)

  expect_false("sparse" %in% res$results_df$id)
})

# ── n_filtered_out component value ───────────────────────────────────────────

test_that("n_filtered_out equals the count of subjects removed by the var/n filter", {
  skip_on_cran()
  set.seed(1)
  base  <- make_panel(n_subjects = 3, n_obs = 25)
  short <- rbind(
    data.frame(id = "s1", time = seq_len(10), x = rnorm(10), y = rnorm(10),
               stringsAsFactors = FALSE),
    data.frame(id = "s2", time = seq_len(10), x = rnorm(10), y = rnorm(10),
               stringsAsFactors = FALSE)
  )
  panel <- rbind(base, short)
  res   <- iarimax(panel, y_series = "y", x_series = "x",
                   id_var = "id", timevar = "time", min_n_subject = 20)
  expect_equal(res$case_number_detail$n_original_df, 5L)
  expect_equal(res$case_number_detail$n_filtered_out, 2L)
  expect_equal(res$case_number_detail$n_used_iarimax, 3L)
})


# ── id_var coercion ───────────────────────────────────────────────────────────

test_that("numeric id_var is coerced to character in results_df", {
  skip_on_cran()
  panel     <- make_panel(n_subjects = 2, n_obs = 25)
  panel$id  <- as.integer(panel$id)

  res <- iarimax(dataframe = panel, y_series = "y", x_series = "x",
                 id_var = "id", timevar = "time")

  expect_type(res$results_df$id, "character")
})

# ── Single valid subject ──────────────────────────────────────────────────────

test_that("single valid subject triggers informative error", {
  skip_on_cran()
  panel <- make_panel(n_subjects = 1, n_obs = 25)

  expect_error(
    iarimax(dataframe = panel, y_series = "y", x_series = "x",
            id_var = "id", timevar = "time"),
    regexp = "not enough cases"
  )
})

# ── Multiple predictors ───────────────────────────────────────────────────────

test_that("two predictors: coefficient columns for both appear in results_df", {
  skip_on_cran()
  panel     <- make_panel(n_subjects = 2, n_obs = 25)
  set.seed(7)
  panel$x2  <- rnorm(nrow(panel))

  res <- iarimax(dataframe = panel, y_series = "y", x_series = c("x", "x2"),
                 focal_predictor = "x", id_var = "id", timevar = "time")

  expect_true("estimate_x"  %in% colnames(res$results_df))
  expect_true("estimate_x2" %in% colnames(res$results_df))
  expect_equal(attr(res, "focal_predictor"), "x")
})

test_that("meta-analysis yi values equal the focal predictor estimates", {
  skip_on_cran()
  panel     <- make_panel(n_subjects = 3, n_obs = 25)
  set.seed(7)
  panel$x2  <- rnorm(nrow(panel))

  res <- iarimax(dataframe = panel, y_series = "y", x_series = c("x", "x2"),
                 focal_predictor = "x", id_var = "id", timevar = "time")

  expect_equal(
    as.numeric(res$meta_analysis$yi),
    res$results_df$estimate_x,
    tolerance = 1e-8
  )
})

# ── Non-sequential timevar ────────────────────────────────────────────────────

test_that("non-sequential timevar (gaps) does not change results vs sequential", {
  skip_on_cran()
  panel     <- make_panel(n_subjects = 2, n_obs = 25)
  panel_gap <- panel
  panel_gap$time <- panel_gap$time * 10

  res_seq <- iarimax(dataframe = panel,     y_series = "y", x_series = "x",
                     id_var = "id", timevar = "time")
  res_gap <- iarimax(dataframe = panel_gap, y_series = "y", x_series = "x",
                     id_var = "id", timevar = "time")

  expect_equal(res_seq$results_df$estimate_x,
               res_gap$results_df$estimate_x, tolerance = 1e-8)
})

# ── Multi-predictor: REMA uses the correct focal predictor column ──────────────

test_that("REMA uses estimate_x2 when focal_predictor is x2", {
  skip_on_cran()
  set.seed(7)
  panel2     <- make_panel(n_subjects = 3, n_obs = 25)
  panel2$x2  <- rnorm(nrow(panel2))

  res <- iarimax(dataframe = panel2, y_series = "y",
                 x_series = c("x", "x2"), focal_predictor = "x2",
                 id_var = "id", timevar = "time")

  valid <- !is.na(res$results_df$estimate_x2)

  expect_equal(
    as.numeric(res$meta_analysis$yi),
    res$results_df$estimate_x2[valid],
    tolerance = 1e-8
  )

  expect_false(isTRUE(all.equal(
    as.numeric(res$meta_analysis$yi),
    res$results_df$estimate_x[valid],
    tolerance = 1e-4
  )))
})

# ── NA in y: Kalman filter produces a valid estimate ─────────────────────────

test_that("subject with some NA in y is still estimated via Kalman filter", {
  skip_on_cran()
  panel_na       <- make_panel(n_subjects = 3, n_obs = 40, seed = 42)
  rows_subj1     <- which(panel_na$id == "1")
  panel_na$y[rows_subj1[1:5]] <- NA

  res  <- iarimax(dataframe = panel_na, y_series = "y", x_series = "x",
                  id_var = "id", timevar = "time")
  row1 <- res$results_df[res$results_df$id == "1", ]

  expect_true("1" %in% res$results_df$id)
  expect_false(is.na(row1$estimate_x))
  expect_true(row1$n_valid < 40)
})


# ── auto.arima failure path ──────────────────────────────────────────────────
# Exercises the tryCatch at the auto.arima call: subject "collinear" has
# x2 = x (perfectly collinear predictors).  Both columns individually pass
# the variance filter but the rank-deficient xreg matrix causes auto.arima
# to error.  The test verifies bookkeeping: the subject appears in
# error_arimax_skipped with NA model stats, but retains a valid raw_cor
# (computed before auto.arima), and the accounting identity holds.

test_that("auto.arima failure: subject in error_arimax_skipped with NA stats and valid raw_cor", {
  skip_on_cran()
  base <- make_panel(n_subjects = 3, n_obs = 25)
  base$x2 <- rnorm(nrow(base))

  set.seed(99)
  bad_x <- rnorm(25)
  bad <- data.frame(id = "collinear", time = seq_len(25),
                    x = bad_x, y = rnorm(25), x2 = bad_x,
                    stringsAsFactors = FALSE)
  panel <- rbind(base, bad)

  res <- suppressMessages(
    iarimax(panel, y_series = "y", x_series = c("x", "x2"),
            focal_predictor = "x", id_var = "id", timevar = "time")
  )

  # Subject passed the filter, so it appears in results_df
  expect_true("collinear" %in% res$results_df$id)

  # Subject should be in error_arimax_skipped
  expect_true("collinear" %in% res$case_number_detail$error_arimax_skipped)

  # Model stats should all be NA
  row <- res$results_df[res$results_df$id == "collinear", ]
  expect_true(is.na(row$nAR))
  expect_true(is.na(row$nI))
  expect_true(is.na(row$nMA))
  expect_true(is.na(row$n_valid))
  expect_true(is.na(row$n_params))
  expect_true(is.na(row$estimate_x))
  expect_true(is.na(row$std.error_x))

  # raw_cor is computed before auto.arima and should survive the failure
  expect_false(is.na(row$raw_cor))

  # Accounting identity must hold
  cd <- res$case_number_detail
  expect_equal(
    cd$n_original_df,
    cd$n_filtered_out + length(cd$error_arimax_skipped) + cd$n_used_iarimax
  )

  # n_used_iarimax excludes the failed subject
  expect_equal(cd$n_used_iarimax, nrow(res$results_df) - length(cd$error_arimax_skipped))
})


# ── Inf values in y or x ────────────────────────────────────────────────────
# Inf is rejected upfront with an informative error.
# NaN (treated as NA by is.na) is handled normally by the variance filter.

test_that("Inf in y triggers an informative error", {
  base <- make_panel(n_subjects = 3, n_obs = 25)
  inf_subj <- data.frame(id = "inf_y", time = seq_len(25),
                         x = rnorm(25), y = c(rnorm(24), Inf),
                         stringsAsFactors = FALSE)
  panel <- rbind(base, inf_subj)

  expect_error(
    iarimax(panel, y_series = "y", x_series = "x",
            id_var = "id", timevar = "time"),
    regexp = "Inf"
  )
})

test_that("Inf in x triggers an informative error", {
  base <- make_panel(n_subjects = 3, n_obs = 25)
  inf_subj <- data.frame(id = "inf_x", time = seq_len(25),
                         x = c(rnorm(24), Inf), y = rnorm(25),
                         stringsAsFactors = FALSE)
  panel <- rbind(base, inf_subj)

  expect_error(
    iarimax(panel, y_series = "y", x_series = "x",
            id_var = "id", timevar = "time"),
    regexp = "Inf"
  )
})

test_that("-Inf in y triggers an informative error", {
  base <- make_panel(n_subjects = 3, n_obs = 25)
  inf_subj <- data.frame(id = "neg_inf", time = seq_len(25),
                         x = rnorm(25), y = c(-Inf, rnorm(24)),
                         stringsAsFactors = FALSE)
  panel <- rbind(base, inf_subj)

  expect_error(
    iarimax(panel, y_series = "y", x_series = "x",
            id_var = "id", timevar = "time"),
    regexp = "Inf"
  )
})

test_that("NaN in y (not NA): treated as missing, filtered if too few remain", {
  skip_on_cran()
  base <- make_panel(n_subjects = 3, n_obs = 25)
  # 15 NaN values → only 10 complete → below min_n_subject = 20
  nan_subj <- data.frame(id = "nan_y", time = seq_len(25),
                         x = rnorm(25), y = c(rep(NaN, 15), rnorm(10)),
                         stringsAsFactors = FALSE)
  panel <- rbind(base, nan_subj)

  res <- iarimax(panel, y_series = "y", x_series = "x",
                 id_var = "id", timevar = "time")

  expect_false("nan_y" %in% res$results_df$id)
})
