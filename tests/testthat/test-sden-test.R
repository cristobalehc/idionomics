# Tests for sden_test()
# Layer 1: fast tests using fake iarimax_results (no auto.arima)
# Layer 2: integration tests on real iarimax output (skip_on_cran)

# ── Helpers ───────────────────────────────────────────────────────────────────

# Fake RMA object: only the fields sden_test reads are needed.
fake_rma <- function(beta, pval) {
  list(beta = matrix(beta, dimnames = list("intrcpt", NULL)), pval = pval)
}

# Fake iarimax_results with known significance counts.
#
# Design (per subject):
#   significant positive  : estimate =  2.0, se = 0.2  -> |t| = 10, df = 48 -> p ≈ 0
#   significant negative  : estimate = -2.0, se = 0.2  -> |t| = 10, df = 48 -> p ≈ 0
#   non-significant       : estimate =  0.1, se = 0.2  -> |t| = 0.5, df = 48 -> p ≈ 0.62
#   failed (NA)           : estimate =  NA,  se = NA   -> pval = NA
#
make_fake_sden_input <- function(n_pos  = 3,
                                 n_neg  = 1,
                                 n_ns   = 2,
                                 n_fail = 0,
                                 rema_beta = 0.5,
                                 rema_pval = 0.03) {
  estimates <- c(rep( 2.0, n_pos),
                 rep(-2.0, n_neg),
                 rep( 0.1, n_ns),
                 rep(  NA, n_fail))
  ses       <- c(rep(0.2, n_pos + n_neg + n_ns),
                 rep( NA, n_fail))
  n_total   <- n_pos + n_neg + n_ns + n_fail

  results_df <- data.frame(
    id            = as.character(seq_len(n_total)),
    estimate_x    = estimates,
    "std.error_x" = ses,
    n_valid       = rep(50L, n_total),
    n_params      = rep(2L,  n_total),
    stringsAsFactors = FALSE,
    check.names   = FALSE
  )

  obj <- list(
    results_df           = results_df,
    meta_analysis        = fake_rma(rema_beta, rema_pval),
    error_arimax_skipped = character(0),
    models               = NULL
  )
  class(obj) <- c("iarimax_results", "list")
  attr(obj, "focal_predictor") <- "x"
  attr(obj, "id_var")          <- "id"
  attr(obj, "timevar")         <- "time"
  obj
}


# ══════════════════════════════════════════════════════════════════════════════
# Layer 1 — Fast tests (no auto.arima, no CRAN skip needed)
# ══════════════════════════════════════════════════════════════════════════════

# ── Validation ────────────────────────────────────────────────────────────────

test_that("sden_test errors on non-iarimax_results input", {
  expect_error(sden_test(list()), "iarimax_results")
})

test_that("sden_test errors on invalid test type", {
  fake <- make_fake_sden_input()
  expect_error(suppressMessages(sden_test(fake, test = "foo")),
               "Invalid test type")
})

test_that("sden_test errors when all p-values are NA", {
  fake <- make_fake_sden_input(n_pos = 0, n_neg = 0, n_ns = 0, n_fail = 3)
  expect_error(suppressMessages(sden_test(fake)), "No valid p-values")
})

test_that("sden_test errors when feature column is missing", {
  fake <- make_fake_sden_input()
  expect_error(suppressMessages(sden_test(fake, feature = "nonexistent")),
               "not found")
})

# ── Output structure ──────────────────────────────────────────────────────────

test_that("sden_test returns sden_results class", {
  fake   <- make_fake_sden_input()
  result <- suppressMessages(sden_test(fake))
  expect_s3_class(result, "sden_results")
})

test_that("sden_test result contains sden_parameters and binomial_test", {
  fake   <- make_fake_sden_input()
  result <- suppressMessages(sden_test(fake))
  expect_true("sden_parameters" %in% names(result))
  expect_true("binomial_test"   %in% names(result))
})

test_that("sden_parameters contains all expected fields", {
  fake     <- make_fake_sden_input()
  result   <- suppressMessages(sden_test(fake))
  expected <- c("test_type", "selection_mechanism", "rema_beta", "pnull",
                "rema_pval", "all_sig_sum", "positive_sig_sum",
                "negative_sig_sum", "number_of_effects", "test_pval",
                "sig_effects")
  expect_true(all(expected %in% names(result$sden_parameters)))
})

test_that("binomial_test is an htest object", {
  fake   <- make_fake_sden_input()
  result <- suppressMessages(sden_test(fake))
  expect_s3_class(result$binomial_test, "htest")
})

test_that("test_pval matches binomial_test p.value", {
  fake   <- make_fake_sden_input()
  result <- suppressMessages(sden_test(fake))
  expect_equal(result$sden_parameters$test_pval,
               result$binomial_test$p.value,
               tolerance = 1e-12)
})

# ── Significance counts ───────────────────────────────────────────────────────

test_that("positive_sig_sum counts subjects with significant positive estimates", {
  fake   <- make_fake_sden_input(n_pos = 3, n_neg = 1, n_ns = 2)
  result <- suppressMessages(sden_test(fake))
  expect_equal(result$sden_parameters$positive_sig_sum, 3L)
})

test_that("negative_sig_sum counts subjects with significant negative estimates", {
  fake   <- make_fake_sden_input(n_pos = 3, n_neg = 1, n_ns = 2)
  result <- suppressMessages(sden_test(fake))
  expect_equal(result$sden_parameters$negative_sig_sum, 1L)
})

test_that("all_sig_sum equals positive_sig_sum + negative_sig_sum", {
  fake   <- make_fake_sden_input(n_pos = 3, n_neg = 1, n_ns = 2)
  result <- suppressMessages(sden_test(fake))
  p      <- result$sden_parameters
  expect_equal(p$all_sig_sum, p$positive_sig_sum + p$negative_sig_sum)
})

test_that("number_of_effects excludes failed subjects (NA estimates)", {
  fake   <- make_fake_sden_input(n_pos = 3, n_neg = 1, n_ns = 2, n_fail = 2)
  result <- suppressMessages(sden_test(fake))
  expect_equal(result$sden_parameters$number_of_effects, 6L)  # 3+1+2, not 8
})

# ── Auto test selection ───────────────────────────────────────────────────────

test_that("auto: non-significant REMA selects ENT", {
  fake   <- make_fake_sden_input(rema_beta = 0.2, rema_pval = 0.40)
  result <- suppressMessages(sden_test(fake, test = "auto"))
  expect_equal(result$sden_parameters$test_type, "ENT")
})

test_that("auto: significant positive REMA selects SDT counter-positive", {
  fake   <- make_fake_sden_input(rema_beta = 0.5, rema_pval = 0.02)
  result <- suppressMessages(sden_test(fake, test = "auto"))
  expect_equal(result$sden_parameters$test_type, "SDT counter-positive")
})

test_that("auto: significant negative REMA selects SDT counter-negative", {
  fake   <- make_fake_sden_input(rema_beta = -0.5, rema_pval = 0.02)
  result <- suppressMessages(sden_test(fake, test = "auto"))
  expect_equal(result$sden_parameters$test_type, "SDT counter-negative")
})

# ── Manual test selection ─────────────────────────────────────────────────────

test_that("test = 'ENT' forces ENT regardless of REMA", {
  # Even with a significant positive REMA, ENT should be forced
  fake   <- make_fake_sden_input(rema_beta = 0.5, rema_pval = 0.02)
  result <- suppressMessages(sden_test(fake, test = "ENT"))
  expect_equal(result$sden_parameters$test_type, "ENT")
  expect_equal(result$sden_parameters$selection_mechanism, "ENT")
})

test_that("test = 'SDT' with positive REMA selects counter-positive", {
  fake   <- make_fake_sden_input(rema_beta = 0.5, rema_pval = 0.40)
  result <- suppressMessages(sden_test(fake, test = "SDT"))
  expect_equal(result$sden_parameters$test_type, "SDT counter-positive")
})

test_that("test = 'SDT' with negative REMA selects counter-negative", {
  fake   <- make_fake_sden_input(rema_beta = -0.5, rema_pval = 0.40)
  result <- suppressMessages(sden_test(fake, test = "SDT"))
  expect_equal(result$sden_parameters$test_type, "SDT counter-negative")
})

# ── pnull correctness ─────────────────────────────────────────────────────────

test_that("ENT uses pnull = alpha_arimax", {
  fake   <- make_fake_sden_input(rema_pval = 0.40)
  result <- suppressMessages(sden_test(fake, alpha_arimax = 0.05))
  expect_equal(result$sden_parameters$pnull, 0.05)
})

test_that("SDT uses pnull = alpha_arimax / 2", {
  fake   <- make_fake_sden_input(rema_beta = 0.5, rema_pval = 0.02)
  result <- suppressMessages(sden_test(fake, alpha_arimax = 0.05))
  expect_equal(result$sden_parameters$pnull, 0.025)
})

test_that("alpha_binom overrides default pnull for ENT", {
  fake    <- make_fake_sden_input(rema_pval = 0.40)
  result  <- suppressMessages(sden_test(fake, alpha_arimax = 0.05,
                                        alpha_binom = 0.10))
  expect_equal(result$sden_parameters$pnull, 0.10)
})

test_that("alpha_binom overrides default pnull for SDT", {
  fake   <- make_fake_sden_input(rema_beta = 0.5, rema_pval = 0.02)
  result <- suppressMessages(sden_test(fake, alpha_arimax = 0.05,
                                       alpha_binom = 0.10))
  expect_equal(result$sden_parameters$pnull, 0.10)
})

test_that("SDT with rema_beta == 0 falls back to ENT behavior", {
  fake   <- make_fake_sden_input(rema_beta = 0, rema_pval = 0.01)
  result <- suppressMessages(sden_test(fake, test = "SDT"))
  # Counts all significant effects (ENT behavior), not just one direction
  expect_equal(result$sden_parameters$test_type, "ENT")
  expect_equal(result$sden_parameters$sig_effects, result$sden_parameters$all_sig_sum)
})

# ── Binomial test correctness ─────────────────────────────────────────────────

test_that("ENT binomial test uses all_sig_sum and alpha_arimax", {
  fake   <- make_fake_sden_input(n_pos = 3, n_neg = 1, n_ns = 2,
                                 rema_pval = 0.40)
  result <- suppressMessages(sden_test(fake, alpha_arimax = 0.05))

  expected <- stats::binom.test(x = 4L, n = 6L, p = 0.05,
                                alternative = "greater")$p.value
  expect_equal(result$sden_parameters$test_pval, expected, tolerance = 1e-12)
})

test_that("SDT counter-positive binomial test uses negative_sig_sum and alpha_arimax/2", {
  fake   <- make_fake_sden_input(n_pos = 3, n_neg = 1, n_ns = 2,
                                 rema_beta = 0.5, rema_pval = 0.02)
  result <- suppressMessages(sden_test(fake, alpha_arimax = 0.05))

  expected <- stats::binom.test(x = 1L, n = 6L, p = 0.025,
                                alternative = "greater")$p.value
  expect_equal(result$sden_parameters$test_pval, expected, tolerance = 1e-12)
})

test_that("SDT counter-negative binomial test uses positive_sig_sum and alpha_arimax/2", {
  fake   <- make_fake_sden_input(n_pos = 3, n_neg = 1, n_ns = 2,
                                 rema_beta = -0.5, rema_pval = 0.02)
  result <- suppressMessages(sden_test(fake, alpha_arimax = 0.05))

  expected <- stats::binom.test(x = 3L, n = 6L, p = 0.025,
                                alternative = "greater")$p.value
  expect_equal(result$sden_parameters$test_pval, expected, tolerance = 1e-12)
})


# ══════════════════════════════════════════════════════════════════════════════
# Layer 2 — Integration on real iarimax output (skip on CRAN)
# ══════════════════════════════════════════════════════════════════════════════

skip_on_cran()

panel  <- make_panel(n_subjects = 4, n_obs = 25, seed = 42)
result <- iarimax(dataframe = panel, y_series = "y", x_series = "x",
                  id_var = "id", timevar = "time")

test_that("sden_test runs without error on real iarimax output", {
  expect_no_error(suppressMessages(sden_test(result)))
})

test_that("sden_test result is sden_results on real output", {
  r <- suppressMessages(sden_test(result))
  expect_s3_class(r, "sden_results")
})

test_that("test_pval matches binomial_test$p.value on real output", {
  r <- suppressMessages(sden_test(result))
  expect_equal(r$sden_parameters$test_pval,
               r$binomial_test$p.value,
               tolerance = 1e-12)
})

test_that("number_of_effects equals subjects with valid estimates on real output", {
  r        <- suppressMessages(sden_test(result))
  n_valid  <- sum(!is.na(result$results_df$estimate_x))
  expect_equal(r$sden_parameters$number_of_effects, n_valid)
})

test_that("sden_test works for non-focal predictor on real multi-predictor output", {
  set.seed(7)
  panel2     <- panel
  panel2$x2  <- rnorm(nrow(panel2))
  result2    <- iarimax(dataframe = panel2, y_series = "y",
                        x_series = c("x", "x2"), focal_predictor = "x",
                        id_var = "id", timevar = "time")

  # Run sden_test on the non-focal predictor x2
  r <- suppressMessages(sden_test(result2, feature = "x2"))
  expect_s3_class(r, "sden_results")
  expect_equal(r$sden_parameters$test_pval,
               r$binomial_test$p.value, tolerance = 1e-12)
})
