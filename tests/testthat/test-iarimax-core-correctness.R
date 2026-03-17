# Verifies that iarimax produces the same numbers as running each
# step manually. These tests are slow (auto.arima runs twice per subject:
# once inside iarimax and once in the manual comparison below) and are
# skipped on CRAN.

skip_on_cran()

panel  <- make_panel(n_subjects = 3, n_obs = 25, seed = 42)
result <- iarimax(dataframe = panel, y_series = "y", x_series = "x",
                  id_var = "id", timevar = "time")

subjects <- unique(panel$id)

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
  for (subj in subjects) {
    m       <- manual_fit(panel, subj)
    pkg_row <- result$results_df[result$results_df$id == subj, ]
    expect_equal(unname(pkg_row$raw_cor), m$cor_p, tolerance = 1e-10,
                 label = paste("raw_cor for subject", subj))
  }
})

test_that("raw_cor (spearman) matches manual cor.test for each subject", {
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
  for (subj in subjects) {
    m       <- manual_fit(panel, subj)
    pkg_row <- result$results_df[result$results_df$id == subj, ]
    expect_equal(unname(pkg_row$n_valid), stats::nobs(m$model),
                 label = paste("n_valid for subject", subj))
  }
})

# ── Temporal ordering ─────────────────────────────────────────────────────────

test_that("shuffled input gives identical results to sorted input", {
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
  df    <- result$results_df
  valid <- !is.na(df[["std.error_x"]])  # metafor silently drops NA rows, so lengths must match
  expect_equal(
    as.numeric(result$meta_analysis$vi),
    df[["std.error_x"]][valid]^2,
    tolerance = 1e-8
  )
})

# ── Signal recovery ───────────────────────────────────────────────────────────
# The true xreg coefficient in make_panel is 0.5. With enough subjects and
# observations the REMA pooled estimate should recover it.

test_that("REMA pooled estimate recovers the true effect (0.5) in a larger panel", {
  big_panel  <- make_panel(n_subjects = 10, n_obs = 35, seed = 42)
  result_big <- iarimax(dataframe = big_panel, y_series = "y", x_series = "x",
                        id_var = "id", timevar = "time")
  rema_est <- as.numeric(result_big$meta_analysis$beta)

  # Correct sign
  expect_true(rema_est > 0,
              label = paste("REMA estimate", round(rema_est, 3), "should be positive"))
  # Within reasonable range of true value 0.5
  expect_true(abs(rema_est - 0.5) < 0.35,
              label = paste("REMA estimate", round(rema_est, 3),
                            "should be within 0.35 of true value 0.5"))
})
