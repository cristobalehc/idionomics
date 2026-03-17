# Edge cases and boundary conditions.
# All tests here run auto.arima, so the whole file is skipped on CRAN.

skip_on_cran()

# ── Filtering: min_n_subject ──────────────────────────────────────────────────

test_that("subject below min_n_subject is absent from results_df", {
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

# ── id_var coercion ───────────────────────────────────────────────────────────

test_that("numeric id_var is coerced to character in results_df", {
  panel     <- make_panel(n_subjects = 2, n_obs = 25)
  panel$id  <- as.integer(panel$id)

  res <- iarimax(dataframe = panel, y_series = "y", x_series = "x",
                 id_var = "id", timevar = "time")

  expect_type(res$results_df$id, "character")
})

# ── Single valid subject ──────────────────────────────────────────────────────

test_that("single valid subject triggers informative error", {
  panel <- make_panel(n_subjects = 1, n_obs = 25)

  expect_error(
    iarimax(dataframe = panel, y_series = "y", x_series = "x",
            id_var = "id", timevar = "time"),
    regexp = "not enough cases"
  )
})

# ── Multiple predictors ───────────────────────────────────────────────────────

test_that("two predictors: coefficient columns for both appear in results_df", {
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
  panel     <- make_panel(n_subjects = 2, n_obs = 25)
  panel_gap <- panel
  panel_gap$time <- panel_gap$time * 10  # e.g. 10, 20, 30, ...

  res_seq <- iarimax(dataframe = panel,     y_series = "y", x_series = "x",
                     id_var = "id", timevar = "time")
  res_gap <- iarimax(dataframe = panel_gap, y_series = "y", x_series = "x",
                     id_var = "id", timevar = "time")

  # Ordering is the same either way; only the scale of time differs
  expect_equal(res_seq$results_df$estimate_x,
               res_gap$results_df$estimate_x, tolerance = 1e-8)
})

# ── Multi-predictor: REMA uses the correct focal predictor column ──────────────

test_that("REMA uses estimate_x2 when focal_predictor is x2", {
  set.seed(7)
  panel2     <- make_panel(n_subjects = 3, n_obs = 25)
  panel2$x2  <- rnorm(nrow(panel2))

  res <- iarimax(dataframe = panel2, y_series = "y",
                 x_series = c("x", "x2"), focal_predictor = "x2",
                 id_var = "id", timevar = "time")

  valid <- !is.na(res$results_df$estimate_x2)

  # REMA yi should match estimate_x2
  expect_equal(
    as.numeric(res$meta_analysis$yi),
    res$results_df$estimate_x2[valid],
    tolerance = 1e-8
  )

  # And should NOT match estimate_x (x and x2 are independent, so estimates differ)
  expect_false(isTRUE(all.equal(
    as.numeric(res$meta_analysis$yi),
    res$results_df$estimate_x[valid],
    tolerance = 1e-4
  )))
})

# ── NA in y: Kalman filter produces a valid estimate ─────────────────────────

test_that("subject with some NA in y is still estimated via Kalman filter", {
  panel_na       <- make_panel(n_subjects = 3, n_obs = 40, seed = 42)
  rows_subj1     <- which(panel_na$id == "1")
  panel_na$y[rows_subj1[1:5]] <- NA  # 35 complete cases remain, above min_n_subject

  res  <- iarimax(dataframe = panel_na, y_series = "y", x_series = "x",
                  id_var = "id", timevar = "time")
  row1 <- res$results_df[res$results_df$id == "1", ]

  # Subject still present and has a valid estimate
  expect_true("1" %in% res$results_df$id)
  expect_false(is.na(row1$estimate_x))

  # nobs reflects the effective observations used (less than 40)
  expect_true(row1$n_valid < 40)
})
