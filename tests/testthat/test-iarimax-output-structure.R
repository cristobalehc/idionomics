# Tests for return value shape, class, and attributes.
# Runs iarimax once (2 subjects, 25 obs) — skipped on CRAN.

skip_on_cran()

panel  <- make_panel(n_subjects = 2, n_obs = 25)
result <- iarimax(dataframe = panel, y_series = "y", x_series = "x",
                  id_var = "id", timevar = "time")

# ── Class and attributes ──────────────────────────────────────────────────────

test_that("result has iarimax_results class", {
  expect_s3_class(result, "iarimax_results")
  expect_true(is.list(result))
})

test_that("focal_predictor attribute matches argument", {
  expect_equal(attr(result, "focal_predictor"), "x")
})

test_that("id_var attribute matches argument", {
  expect_equal(attr(result, "id_var"), "id")
})

test_that("timevar attribute matches argument", {
  expect_equal(attr(result, "timevar"), "time")
})

# ── List structure ────────────────────────────────────────────────────────────

test_that("result has exactly the expected top-level names", {
  expect_setequal(
    names(result),
    c("results_df", "meta_analysis", "case_number_detail", "models")
  )
})

test_that("case_number_detail contains all expected fields", {
  expected <- c("n_original_df", "n_used_iarimax", "n_filtered_out", "error_arimax_skipped")
  expect_true(all(expected %in% names(result$case_number_detail)))
})

test_that("models is NULL when keep_models = FALSE", {
  expect_null(result$models)
})

test_that("models is a named list of Arima objects when keep_models = TRUE", {
  res <- iarimax(dataframe = panel, y_series = "y", x_series = "x",
                 id_var = "id", timevar = "time", keep_models = TRUE)
  expect_type(res$models, "list")
  expect_length(res$models, nrow(res$results_df))
  valid_models <- Filter(Negate(is.null), res$models)
  # "ARIMA" (fable class) is stripped; "Arima" (forecast class) is kept
  expect_true(all(sapply(valid_models, function(m) "Arima" %in% class(m))))
})

test_that("error_arimax_skipped is a character vector inside case_number_detail", {
  expect_type(result$case_number_detail$error_arimax_skipped, "character")
})

test_that("meta_analysis is a metafor rma object", {
  expect_s3_class(result$meta_analysis, "rma")
})

# ── results_df columns ────────────────────────────────────────────────────────

test_that("results_df contains all required columns", {
  required <- c("id", "nAR", "nI", "nMA", "raw_cor",
                "n_valid", "n_params", "estimate_x", "std.error_x")
  expect_true(all(required %in% colnames(result$results_df)))
})

test_that("results_df id column is character type", {
  expect_type(result$results_df$id, "character")
})

test_that("results_df has one row per subject that passed filtering", {
  expect_equal(nrow(result$results_df), length(unique(panel$id)))
})

test_that("n_valid is non-negative for all fitted subjects", {
  valid_rows <- result$results_df[!is.na(result$results_df$n_valid), ]
  expect_true(all(valid_rows$n_valid >= 0))
})

test_that("verbose = TRUE produces same numeric results as verbose = FALSE", {
  res_v <- suppressMessages(
    iarimax(dataframe = panel, y_series = "y", x_series = "x",
            id_var = "id", timevar = "time", verbose = TRUE)
  )
  expect_equal(result$results_df$estimate_x,  res_v$results_df$estimate_x)
  expect_equal(result$results_df$std.error_x, res_v$results_df$std.error_x)
  expect_equal(result$results_df$raw_cor,     res_v$results_df$raw_cor)
})
