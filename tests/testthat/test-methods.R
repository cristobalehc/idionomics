# Tests for S3 methods and downstream behavior of NULL meta_analysis.
# All Layer 1 (fake objects, no auto.arima, no skip_on_cran needed).

# ── Helpers ───────────────────────────────────────────────────────────────────

make_fake_sden_result <- function(test_type = "ENT",
                                  selection_mechanism = "auto",
                                  rema_beta = 0.2,
                                  rema_pval = 0.40,
                                  all_sig = 4L, pos = 3L, neg = 1L,
                                  n_effects = 6L, sig_effects = 4L,
                                  pnull = 0.05) {
  btest <- stats::binom.test(x = sig_effects, n = n_effects,
                             p = pnull, alternative = "greater")

  params <- list(
    test_type          = test_type,
    selection_mechanism = selection_mechanism,
    rema_beta          = rema_beta,
    pnull              = pnull,
    rema_pval          = rema_pval,
    all_sig_sum        = all_sig,
    positive_sig_sum   = pos,
    negative_sig_sum   = neg,
    number_of_effects  = n_effects,
    test_pval          = btest$p.value,
    sig_effects        = sig_effects
  )

  result <- list(sden_parameters = params, binomial_test = btest)
  class(result) <- c("sden_results", "list")
  result
}

# ── summary.sden_results: ENT / auto ─────────────────────────────────────────

test_that("summary.sden_results runs without error for ENT (auto)", {
  r <- make_fake_sden_result()
  expect_no_error(capture.output(summary(r)))
})

test_that("summary.sden_results output contains test type label for ENT", {
  r <- make_fake_sden_result()
  expect_output(summary(r), regexp = "ENT")
})

test_that("summary.sden_results output contains 'Automatic' for auto selection", {
  r <- make_fake_sden_result(selection_mechanism = "auto")
  expect_output(summary(r), regexp = "Automatic")
})

test_that("summary.sden_results output contains REMA beta", {
  r <- make_fake_sden_result(rema_beta = 0.2)
  expect_output(summary(r), regexp = "0.2")
})

test_that("summary.sden_results output contains the binomial p-value", {
  r <- make_fake_sden_result()
  expect_output(summary(r), regexp = "p-value")
})

test_that("summary.sden_results returns the object invisibly", {
  r <- make_fake_sden_result()
  captured <- NULL
  capture.output(captured <- withVisible(summary(r)))
  expect_false(captured$visible)
  expect_s3_class(captured$value, "sden_results")
})

# ── summary.sden_results: SDT counter-positive / manual ──────────────────────

test_that("summary.sden_results output contains 'SDT' for SDT counter-positive", {
  r <- make_fake_sden_result(test_type = "SDT counter-positive",
                             selection_mechanism = "SDT",
                             rema_beta = 0.5, rema_pval = 0.02,
                             sig_effects = 1L, pnull = 0.025)
  expect_output(summary(r), regexp = "SDT")
})

test_that("summary.sden_results output contains 'Manual' for manual selection", {
  r <- make_fake_sden_result(test_type = "SDT counter-positive",
                             selection_mechanism = "SDT",
                             rema_beta = 0.5, rema_pval = 0.02,
                             sig_effects = 1L, pnull = 0.025)
  expect_output(summary(r), regexp = "Manual")
})

# ── Helpers for iarimax_results S3 methods ────────────────────────────────────

make_fake_iarimax <- function(with_meta = TRUE) {
  df <- data.frame(
    id            = c("1", "2", "3", "4"),
    nAR           = c(1L, 0L, 1L, 0L),
    nI            = c(0L, 0L, 1L, 0L),
    nMA           = c(0L, 1L, 0L, 1L),
    raw_cor       = c(0.3, 0.5, -0.1, 0.2),
    n_valid       = c(23L, 24L, 22L, 25L),
    n_params      = c(2L,  2L,  3L,  2L),
    estimate_x    = c(0.4,  0.6, -0.2,  0.3),
    `std.error_x` = c(0.10, 0.12, 0.08, 0.11),
    stringsAsFactors = FALSE,
    check.names      = FALSE
  )

  ma <- if (with_meta) {
    suppressWarnings(metafor::rma(yi = df$estimate_x, sei = df[["std.error_x"]], method = "REML"))
  } else {
    NULL
  }

  result <- list(
    results_df     = df,
    meta_analysis  = ma,
    case_number_detail = list(
      n_original_df        = 6L,
      n_filtered_out       = 2L,
      error_arimax_skipped = character(0),
      n_used_iarimax       = 4L
    ),
    models = NULL
  )
  class(result) <- c("iarimax_results", "list")
  attr(result, "outcome")         <- "y"
  attr(result, "focal_predictor") <- "x"
  attr(result, "id_var")          <- "id"
  attr(result, "timevar")         <- "time"
  result
}

# ── summary.iarimax_results: Layer 1 (fake object) ───────────────────────────

test_that("summary.iarimax_results runs without error", {
  r <- make_fake_iarimax()
  expect_no_error(capture.output(summary(r)))
})

test_that("summary.iarimax_results runs without error when meta_analysis is NULL", {
  r <- make_fake_iarimax(with_meta = FALSE)
  expect_no_error(capture.output(summary(r)))
})

test_that("summary.iarimax_results returns object invisibly", {
  r <- make_fake_iarimax()
  captured <- NULL
  capture.output(captured <- withVisible(summary(r)))
  expect_false(captured$visible)
  expect_s3_class(captured$value, "iarimax_results")
})

test_that("summary.iarimax_results output contains focal predictor name", {
  r <- make_fake_iarimax()
  expect_output(summary(r), regexp = "x")
})

test_that("summary.iarimax_results output contains n_original count", {
  r <- make_fake_iarimax()
  expect_output(summary(r), regexp = "6")
})

test_that("summary.iarimax_results output contains REMA beta when meta_analysis present", {
  r <- make_fake_iarimax()
  expect_output(summary(r), regexp = "beta")
})

test_that("summary.iarimax_results output contains heterogeneity section when meta_analysis present", {
  r <- make_fake_iarimax()
  expect_output(summary(r), regexp = "tau2")
})

test_that("summary.iarimax_results output contains 'not available' when meta_analysis is NULL", {
  r <- make_fake_iarimax(with_meta = FALSE)
  expect_output(summary(r), regexp = "not available")
})

test_that("summary.iarimax_results respects custom alpha argument", {
  r <- make_fake_iarimax()
  expect_output(summary(r, alpha = 0.01), regexp = "0.01")
})

# ── plot.iarimax_results: guard clause (Layer 1) ──────────────────────────────

test_that("plot.iarimax_results errors on unknown feature naming the bad variable", {
  r <- make_fake_iarimax()
  expect_error(plot(r, feature = "not_a_feature"), regexp = "not_a_feature")
})

# ── summary.iarimax_results + plot.iarimax_results: Layer 2 (real iarimax) ───

test_that("summary.iarimax_results runs without error on real iarimax output", {
  skip_on_cran()
  panel  <- make_panel(n_subjects = 4, n_obs = 25)
  result <- iarimax(dataframe = panel, y_series = "y", x_series = "x",
                    id_var = "id", timevar = "time")
  expect_no_error(capture.output(summary(result)))
})

test_that("summary.iarimax_results reports correct n_original on real output", {
  skip_on_cran()
  panel  <- make_panel(n_subjects = 4, n_obs = 25)
  result <- iarimax(dataframe = panel, y_series = "y", x_series = "x",
                    id_var = "id", timevar = "time")
  expect_output(summary(result), regexp = as.character(length(unique(panel$id))))
})

test_that("plot.iarimax_results returns a ggplot object on real output", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  panel  <- make_panel(n_subjects = 4, n_obs = 25)
  result <- iarimax(dataframe = panel, y_series = "y", x_series = "x",
                    id_var = "id", timevar = "time")
  expect_s3_class(plot(result), "ggplot")
})

test_that("plot.iarimax_results works with custom lims", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  panel  <- make_panel(n_subjects = 4, n_obs = 25)
  result <- iarimax(dataframe = panel, y_series = "y", x_series = "x",
                    id_var = "id", timevar = "time")
  expect_no_error(plot(result, lims = c(-2, 2)))
})

test_that("plot.iarimax_results works with substantive name labels", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  panel  <- make_panel(n_subjects = 4, n_obs = 25)
  result <- iarimax(dataframe = panel, y_series = "y", x_series = "x",
                    id_var = "id", timevar = "time")
  expect_s3_class(plot(result, y_series_name = "Mood", x_series_name = "Stress"), "ggplot")
})

test_that("plot.iarimax_results still produces a plot when meta_analysis is NULL", {
  skip_if_not_installed("ggplot2")
  r <- make_fake_iarimax(with_meta = FALSE)
  expect_s3_class(plot(r), "ggplot")
})

# ── sden_test errors when meta_analysis is NULL ───────────────────────────────

test_that("sden_test stops with informative error when meta_analysis is NULL", {
  # Simulates the case where too few subjects produced valid models and
  # metafor::rma() failed, leaving meta_analysis = NULL in the iarimax object.
  fake <- list(
    results_df = data.frame(
      id            = c("1", "2"),
      estimate_x    = c(0.5, -0.3),
      "std.error_x" = c(0.1,  0.15),
      n_valid       = c(25L,  25L),
      n_params      = c(2L,   2L),
      stringsAsFactors = FALSE,
      check.names   = FALSE
    ),
    meta_analysis        = NULL,
    error_arimax_skipped = character(0),
    models               = NULL
  )
  class(fake) <- c("iarimax_results", "list")
  attr(fake, "focal_predictor") <- "x"
  attr(fake, "id_var")          <- "id"
  attr(fake, "timevar")         <- "time"

  expect_error(suppressMessages(sden_test(fake)), regexp = "SDEN test stopped")
})
