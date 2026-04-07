# Verifies that iarimax produces the same numbers as running each
# step manually. These tests are slow (auto.arima runs twice per subject:
# once inside iarimax and once in the manual comparison below) and are
# skipped on CRAN.
#
# Uses the same memoised accessor pattern as test-iarimax-output-structure.R
# so that skip_on_cran() fires inside each test, not at file level.

panel    <- make_panel(n_subjects = 3, n_obs = 25, seed = 42)
subjects <- unique(panel$id)

.cc_cache <- new.env(parent = emptyenv())

.get_result <- function() {
  if (!exists("r", envir = .cc_cache)) {
    skip_on_cran()
    .cc_cache$r <- iarimax(dataframe = panel, y_series = "y", x_series = "x",
                           id_var = "id", timevar = "time")
  }
  .cc_cache$r
}

# Helper: reproduce exactly what iarimax does for one subject.
manual_fit <- function(data, subj) {
  sub   <- data[data$id == subj, ]
  sub   <- sub[order(sub$time), ]
  y_vec <- sub$y
  x_mat <- as.matrix(sub[, "x", drop = FALSE])

  model <- forecast::auto.arima(y = y_vec, xreg = x_mat,
                                approximation = FALSE, stepwise = FALSE)
  class(model) <- setdiff(class(model), "ARIMA")  # mirrors iarimax's workaround

  list(
    model  = model,
    tidy   = broom::tidy(model),
    cor_p  = suppressWarnings(
      stats::cor.test(sub$x, sub$y, method = "pearson")
    )$estimate[[1]],
    cor_s  = suppressWarnings(
      stats::cor.test(sub$x, sub$y, method = "spearman")
    )$estimate[[1]]
  )
}

# ── xreg coefficient ─────────────────────────────────────────────────────────

test_that("estimate_x matches manual auto.arima for each subject", {
  result <- .get_result()
  for (subj in subjects) {
    m       <- manual_fit(panel, subj)
    pkg_row <- result$results_df[result$results_df$id == subj, ]
    expect_equal(
      pkg_row$estimate_x,
      m$tidy$estimate[m$tidy$term == "x"],
      tolerance = 1e-8,
      label = paste("estimate_x for subject", subj)
    )
  }
})

test_that("std.error_x matches manual auto.arima for each subject", {
  result <- .get_result()
  for (subj in subjects) {
    m       <- manual_fit(panel, subj)
    pkg_row <- result$results_df[result$results_df$id == subj, ]
    expect_equal(
      pkg_row$std.error_x,
      m$tidy$std.error[m$tidy$term == "x"],
      tolerance = 1e-8,
      label = paste("std.error_x for subject", subj)
    )
  }
})

# ── ARIMA orders ──────────────────────────────────────────────────────────────

test_that("ARIMA order (nAR, nI, nMA) matches manual auto.arima for each subject", {
  result <- .get_result()
  for (subj in subjects) {
    m       <- manual_fit(panel, subj)
    pkg_row <- result$results_df[result$results_df$id == subj, ]
    expect_equal(unname(pkg_row$nAR), m$model$arma[1], label = paste("nAR for subject", subj))
    expect_equal(unname(pkg_row$nI),  m$model$arma[6], label = paste("nI  for subject", subj))
    expect_equal(unname(pkg_row$nMA), m$model$arma[2], label = paste("nMA for subject", subj))
  }
})

# ── Raw correlation ───────────────────────────────────────────────────────────

test_that("raw_cor (pearson) matches manual cor.test for each subject", {
  result <- .get_result()
  for (subj in subjects) {
    m       <- manual_fit(panel, subj)
    pkg_row <- result$results_df[result$results_df$id == subj, ]
    expect_equal(unname(pkg_row$raw_cor), m$cor_p, tolerance = 1e-10,
                 label = paste("raw_cor for subject", subj))
  }
})

test_that("raw_cor (spearman) matches manual cor.test for each subject", {
  skip_on_cran()
  result_s <- iarimax(dataframe = panel, y_series = "y", x_series = "x",
                      id_var = "id", timevar = "time",
                      correlation_method = "spearman")
  for (subj in subjects) {
    m       <- manual_fit(panel, subj)
    pkg_row <- result_s$results_df[result_s$results_df$id == subj, ]
    expect_equal(unname(pkg_row$raw_cor), m$cor_s, tolerance = 1e-10,
                 label = paste("spearman raw_cor for subject", subj))
  }
})

test_that("raw_cor (kendall) matches manual cor.test for each subject", {
  skip_on_cran()
  result_k <- iarimax(dataframe = panel, y_series = "y", x_series = "x",
                      id_var = "id", timevar = "time",
                      correlation_method = "kendall")
  for (subj in subjects) {
    expected <- suppressWarnings(
      stats::cor.test(panel$x[panel$id == subj],
                      panel$y[panel$id == subj],
                      method = "kendall")
    )$estimate[[1]]
    pkg_row <- result_k$results_df[result_k$results_df$id == subj, ]
    expect_equal(unname(pkg_row$raw_cor), expected, tolerance = 1e-10,
                 label = paste("kendall raw_cor for subject", subj))
  }
})

# ── n_valid ───────────────────────────────────────────────────────────────────

test_that("n_valid matches stats::nobs of manual model for each subject", {
  result <- .get_result()
  for (subj in subjects) {
    m       <- manual_fit(panel, subj)
    pkg_row <- result$results_df[result$results_df$id == subj, ]
    expect_equal(unname(pkg_row$n_valid), stats::nobs(m$model),
                 label = paste("n_valid for subject", subj))
  }
})

# ── Temporal ordering ─────────────────────────────────────────────────────────

test_that("shuffled input gives identical results to sorted input", {
  skip_on_cran()
  result  <- .get_result()
  set.seed(99)
  shuffled        <- panel[sample(nrow(panel)), ]
  result_shuffled <- iarimax(dataframe = shuffled, y_series = "y", x_series = "x",
                             id_var = "id", timevar = "time")

  r1 <- result$results_df[order(result$results_df$id), ]
  r2 <- result_shuffled$results_df[order(result_shuffled$results_df$id), ]

  expect_equal(r1$estimate_x,  r2$estimate_x,  tolerance = 1e-8)
  expect_equal(r1$std.error_x, r2$std.error_x, tolerance = 1e-8)
  expect_equal(r1$raw_cor,     r2$raw_cor,     tolerance = 1e-10)
  expect_equal(r1$nAR,         r2$nAR)
  expect_equal(r1$nI,          r2$nI)
  expect_equal(r1$nMA,         r2$nMA)
})

# ── n_params ──────────────────────────────────────────────────────────────────

test_that("n_params matches length(model$coef) for each subject", {
  result <- .get_result()
  for (subj in subjects) {
    m       <- manual_fit(panel, subj)
    pkg_row <- result$results_df[result$results_df$id == subj, ]
    expect_equal(
      unname(pkg_row$n_params),
      length(m$model$coef),
      label = paste("n_params for subject", subj)
    )
  }
})

# ── REMA standard errors ──────────────────────────────────────────────────────

test_that("REMA sampling variances equal std.error_x squared", {
  result <- .get_result()
  df     <- result$results_df
  valid  <- !is.na(df[["std.error_x"]])
  expect_equal(
    as.numeric(result$meta_analysis$vi),
    df[["std.error_x"]][valid]^2,
    tolerance = 1e-8
  )
})

# ── fixed_d ───────────────────────────────────────────────────────────────────

test_that("fixed_d = 0 forces nI = 0 for all subjects", {
  skip_on_cran()
  result_d0 <- iarimax(dataframe = panel, y_series = "y", x_series = "x",
                       id_var = "id", timevar = "time", fixed_d = 0)
  expect_true(all(result_d0$results_df$nI == 0, na.rm = TRUE))
})

test_that("fixed_d = 1 forces nI = 1 for all subjects", {
  skip_on_cran()
  result_d1 <- iarimax(dataframe = panel, y_series = "y", x_series = "x",
                       id_var = "id", timevar = "time", fixed_d = 1)
  expect_true(all(result_d1$results_df$nI == 1, na.rm = TRUE))
})

test_that("fixed_d = NULL allows varying nI across subjects", {
  result <- .get_result()
  # This is not guaranteed to vary with 3 subjects, but we verify it runs
  # and nI contains valid values
  expect_true(all(result$results_df$nI >= 0, na.rm = TRUE))
})

# ── Signal recovery ───────────────────────────────────────────────────────────
# The true xreg coefficient in make_panel is 0.5. With enough subjects and
# observations the REMA pooled estimate should recover it.

# ── REMA independent verification ────────────────────────────────────────────
# Verifies that the meta_analysis object stored by iarimax matches an
# independent metafor::rma() call on the same data (beta, CI, tau2, I2).

test_that("REMA beta matches independent metafor::rma() computation", {
  result <- .get_result()
  df     <- result$results_df
  valid  <- !is.na(df[["estimate_x"]])

  manual <- metafor::rma(yi = df[["estimate_x"]][valid],
                         sei = df[["std.error_x"]][valid],
                         method = "REML")

  expect_equal(as.numeric(result$meta_analysis$beta),
               as.numeric(manual$beta), tolerance = 1e-10)
})

test_that("REMA CI bounds match independent metafor::rma() computation", {
  result <- .get_result()
  df     <- result$results_df
  valid  <- !is.na(df[["estimate_x"]])

  manual <- metafor::rma(yi = df[["estimate_x"]][valid],
                         sei = df[["std.error_x"]][valid],
                         method = "REML")

  expect_equal(result$meta_analysis$ci.lb, manual$ci.lb, tolerance = 1e-10)
  expect_equal(result$meta_analysis$ci.ub, manual$ci.ub, tolerance = 1e-10)
})

test_that("REMA tau2 and I2 match independent metafor::rma() computation", {
  result <- .get_result()
  df     <- result$results_df
  valid  <- !is.na(df[["estimate_x"]])

  manual <- metafor::rma(yi = df[["estimate_x"]][valid],
                         sei = df[["std.error_x"]][valid],
                         method = "REML")

  expect_equal(result$meta_analysis$tau2, manual$tau2, tolerance = 1e-10)
  expect_equal(result$meta_analysis$I2,   manual$I2,   tolerance = 1e-8)
})

# ── Signal recovery ───────────────────────────────────────────────────────────

test_that("REMA pooled estimate recovers the true effect (0.5) in a larger panel", {
  skip_on_cran()
  big_panel  <- make_panel(n_subjects = 10, n_obs = 35, seed = 42)
  result_big <- iarimax(dataframe = big_panel, y_series = "y", x_series = "x",
                        id_var = "id", timevar = "time")
  rema_est <- as.numeric(result_big$meta_analysis$beta)

  expect_true(rema_est > 0,
              label = paste("REMA estimate", round(rema_est, 3), "should be positive"))
  expect_true(abs(rema_est - 0.5) < 0.35,
              label = paste("REMA estimate", round(rema_est, 3),
                            "should be within 0.35 of true value 0.5"))
})

# ── REML method guard ─────────────────────────────────────────────────────────

test_that("meta_analysis uses REML estimation method", {
  result <- .get_result()
  expect_match(result$meta_analysis$method, "REML")
})

# ── End-to-end pipeline integration ──────────────────────────────────────────

test_that("full pipeline (i_screener -> pmstandardize -> i_detrender -> iarimax -> i_pval -> sden_test) runs end-to-end", {
  skip_on_cran()
  panel <- make_panel(n_subjects = 5, n_obs = 35, seed = 123)

  # Step 0: screen
  df_clean <- i_screener(panel, cols = c("y", "x"), id_var = "id",
                         min_n_subject = 20)

  # Step 1: standardize
  df_psd <- pmstandardize(df_clean, cols = c("y", "x"), id_var = "id")
  expect_true("y_psd" %in% names(df_psd))
  expect_true("x_psd" %in% names(df_psd))

  # Step 2: detrend
  df_ready <- i_detrender(df_psd, cols = c("y_psd", "x_psd"),
                          id_var = "id", timevar = "time")
  expect_true("y_psd_dt" %in% names(df_ready))
  expect_true("x_psd_dt" %in% names(df_ready))

  # Step 3: model
  result <- iarimax(df_ready, y_series = "y_psd_dt", x_series = "x_psd_dt",
                    id_var = "id", timevar = "time")
  expect_s3_class(result, "iarimax_results")
  expect_true(nrow(result$results_df) >= 2)

  # Step 4: p-values
  result <- i_pval(result)
  expect_true("pval_x_psd_dt" %in% names(result$results_df))

  # Step 5: sden_test
  sden <- suppressMessages(sden_test(result))
  expect_s3_class(sden, "sden_results")
  expect_true(sden$sden_parameters$test_pval >= 0)
  expect_true(sden$sden_parameters$test_pval <= 1)
})
