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

test_that("summary.sden_results auto + SDT counter-positive prints correct explanation", {
  # Exercises lines 228-231 in summary.sden_results: auto selection with a
  # positive-and-significant pooled effect.
  r <- make_fake_sden_result(test_type = "SDT counter-positive",
                             selection_mechanism = "auto",
                             rema_beta = 0.5, rema_pval = 0.02,
                             sig_effects = 1L, pnull = 0.025)
  expect_output(summary(r), regexp = "positive")
  expect_output(summary(r), regexp = "Automatic")
})

test_that("summary.sden_results auto + SDT counter-negative prints correct explanation and hypothesis", {
  # Exercises lines 231-234 (explanation) and lines 248-250 (hypothesis) in
  # summary.sden_results: auto selection with a negative-and-significant effect.
  r <- make_fake_sden_result(test_type = "SDT counter-negative",
                             selection_mechanism = "auto",
                             rema_beta = -0.5, rema_pval = 0.02,
                             sig_effects = 3L, pnull = 0.025)
  expect_output(summary(r), regexp = "negative")
  expect_output(summary(r), regexp = "Automatic")
  expect_output(summary(r), regexp = "positive significant")
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

# ── summary.iarimax_results: direction/significance counts (Layer 1) ────────

test_that("summary.iarimax_results prints correct direction and significance counts", {
  # Fake data: estimates = 0.4, 0.6, -0.2, 0.3; SEs = 0.10, 0.12, 0.08, 0.11
  # df = n_valid - n_params: 21, 22, 19, 23
  # t-stats: 4.0, 5.0, -2.5, 2.73
  # All have |t| > qt(0.975, df) ≈ 2.08, so all 4 are significant at alpha=0.05
  # Directions: 3 positive (subjects 1,2,4), 1 negative (subject 3)
  r   <- make_fake_iarimax()
  out <- capture.output(summary(r))
  combined <- paste(out, collapse = " ")
  expect_match(combined, "Positive\\s*:\\s*3")
  expect_match(combined, "Significant positive\\s*:\\s*3")
  expect_match(combined, "Negative\\s*:\\s*1")
  expect_match(combined, "Significant negative\\s*:\\s*1")
})

test_that("summary.iarimax_results respects alpha for significance counts", {
  # At alpha = 0.001, subject 3 (|t|=2.5, df=19) should no longer be significant
  # qt(0.9995, 19) ≈ 3.88 > 2.5
  # Subjects 1 (|t|=4.0, df=21) and 2 (|t|=5.0, df=22) remain significant
  # Subject 4 (|t|=2.73, df=23): qt(0.9995, 23) ≈ 3.77 > 2.73, not significant
  r   <- make_fake_iarimax()
  out <- capture.output(summary(r, alpha = 0.001))
  combined <- paste(out, collapse = " ")
  # Directions unchanged: still 3 positive, 1 negative
  expect_match(combined, "Positive\\s*:\\s*3")
  expect_match(combined, "Negative\\s*:\\s*1")
  # Sig positive: only subjects 1 and 2
  expect_match(combined, "Significant positive\\s*:\\s*2")
  # Sig negative: subject 3 no longer significant
  expect_match(combined, "Significant negative\\s*:\\s*0")
})

# ── summary.iarimax_results: alpha validation (Layer 1) ─────────────────────

test_that("summary.iarimax_results errors on non-numeric alpha", {
  r <- make_fake_iarimax()
  expect_error(capture.output(summary(r, alpha = "x")), regexp = "alpha")
})

test_that("summary.iarimax_results errors on alpha = 0", {
  r <- make_fake_iarimax()
  expect_error(capture.output(summary(r, alpha = 0)), regexp = "alpha")
})

test_that("summary.iarimax_results errors on alpha = 1", {
  r <- make_fake_iarimax()
  expect_error(capture.output(summary(r, alpha = 1)), regexp = "alpha")
})

test_that("summary.iarimax_results errors on alpha = Inf", {
  r <- make_fake_iarimax()
  expect_error(capture.output(summary(r, alpha = Inf)), regexp = "alpha")
})

# ── plot.iarimax_results: guard clause (Layer 1) ──────────────────────────────

test_that("plot.iarimax_results errors on non-numeric alpha_crit_t", {
  r <- make_fake_iarimax()
  expect_error(plot(r, alpha_crit_t = "x"), regexp = "alpha_crit_t")
})

test_that("plot.iarimax_results errors on alpha_crit_t = 0", {
  r <- make_fake_iarimax()
  expect_error(plot(r, alpha_crit_t = 0), regexp = "alpha_crit_t")
})

test_that("plot.iarimax_results errors on alpha_crit_t = 1", {
  r <- make_fake_iarimax()
  expect_error(plot(r, alpha_crit_t = 1), regexp = "alpha_crit_t")
})

test_that("plot.iarimax_results errors on non-numeric lims", {
  r <- make_fake_iarimax()
  expect_error(plot(r, lims = "bad"), regexp = "lims")
})

test_that("plot.iarimax_results errors on lims of wrong length", {
  r <- make_fake_iarimax()
  expect_error(plot(r, lims = c(-1, 0, 1)), regexp = "lims")
})

test_that("plot.iarimax_results errors on lims with Inf", {
  r <- make_fake_iarimax()
  expect_error(plot(r, lims = c(-Inf, 1)), regexp = "lims")
})

test_that("plot.iarimax_results errors on unknown feature naming the bad variable", {
  r <- make_fake_iarimax()
  expect_error(plot(r, feature = "not_a_feature"), regexp = "not_a_feature")
})

test_that("plot.iarimax_results works for a non-focal feature present in results_df", {
  skip_if_not_installed("ggplot2")
  r <- make_fake_iarimax()
  r$results_df$estimate_z     <- c(0.1, -0.2, 0.3, 0.0)
  r$results_df[["std.error_z"]] <- c(0.05, 0.06, 0.07, 0.04)
  expect_s3_class(plot(r, feature = "z"), "ggplot")
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

test_that("plot.iarimax_results y-axis and caption reflect custom alpha_crit_t", {
  skip_if_not_installed("ggplot2")
  r <- make_fake_iarimax()
  p <- plot(r, alpha_crit_t = 0.01)
  expect_match(p$labels$y,       "99%")
  expect_match(p$labels$caption, "99%")
})

test_that("plot.iarimax_results default alpha_crit_t uses 95% in y-axis and caption", {
  skip_if_not_installed("ggplot2")
  r <- make_fake_iarimax()
  p <- plot(r)
  expect_match(p$labels$y,       "95%")
  expect_match(p$labels$caption, "95%")
})

test_that("plot.iarimax_results line_color maps correctly to significance", {
  skip_if_not_installed("ggplot2")
  r  <- make_fake_iarimax()
  p  <- plot(r)
  pd <- ggplot2::ggplot_build(p)

  # The linerange layer is the second geom (after geom_point)
  line_data <- pd$data[[2]]

  # Subject 3: estimate = -0.2, SE = 0.08, df = 22 - 3 = 19
  # CI entirely negative -> should be red
  # Subject 1: estimate = 0.4, SE = 0.10, df = 23 - 2 = 21
  # CI: 0.4 +/- qt(0.975,21)*0.10 = [0.19, 0.61] -> entirely positive -> green
  expect_true("red"   %in% line_data$colour)
  expect_true("green" %in% line_data$colour)
})

test_that("plot.iarimax_results REMA band matches meta_analysis CI", {
  skip_if_not_installed("ggplot2")
  r  <- make_fake_iarimax()
  p  <- plot(r)
  pd <- ggplot2::ggplot_build(p)

  # The annotate("rect") layer is the last one (after point, linerange, hline, hline)
  rect_data <- pd$data[[length(pd$data)]]
  expect_equal(rect_data$ymin[1], r$meta_analysis$ci.lb, tolerance = 1e-8)
  expect_equal(rect_data$ymax[1], r$meta_analysis$ci.ub, tolerance = 1e-8)
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

  expect_error(suppressMessages(sden_test(fake)), regexp = "meta_analysis")
})
