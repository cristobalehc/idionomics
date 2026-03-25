# Tests for return value shape, class, and attributes.
#
# Why the accessor pattern instead of a file-level skip_on_cran()?
# A file-level skip_on_cran() aborts the entire file on CRAN (and under
# covr), so a single iarimax() failure silently drops all tests with no
# per-test diagnostic. The accessor below runs iarimax() once, caches the
# result, and returns it on every subsequent call. skip_on_cran() lives
# inside the accessor, so it fires only when a test actually calls it —
# each test is skipped individually, and file-level code (panel creation,
# function definitions) always runs cleanly.

panel <- make_panel(n_subjects = 2, n_obs = 25)

.oc_cache <- new.env(parent = emptyenv())

.get_result <- function() {
  if (!exists("r", envir = .oc_cache)) {
    skip_on_cran()
    .oc_cache$r <- iarimax(dataframe = panel, y_series = "y", x_series = "x",
                           id_var = "id", timevar = "time")
  }
  .oc_cache$r
}


# ── Class and attributes ──────────────────────────────────────────────────────

test_that("result has iarimax_results class", {
  result <- .get_result()
  expect_s3_class(result, "iarimax_results")
  expect_true(is.list(result))
})

test_that("focal_predictor attribute matches argument", {
  expect_equal(attr(.get_result(), "focal_predictor"), "x")
})

test_that("id_var attribute matches argument", {
  expect_equal(attr(.get_result(), "id_var"), "id")
})

test_that("timevar attribute matches argument", {
  expect_equal(attr(.get_result(), "timevar"), "time")
})

test_that("outcome attribute matches y_series argument", {
  expect_equal(attr(.get_result(), "outcome"), "y")
})


# ── List structure ────────────────────────────────────────────────────────────

test_that("result has exactly the expected top-level names", {
  expect_setequal(
    names(.get_result()),
    c("results_df", "meta_analysis", "case_number_detail", "models")
  )
})

test_that("case_number_detail contains all expected fields", {
  expected <- c("n_original_df", "n_used_iarimax", "n_filtered_out",
                "error_arimax_skipped")
  expect_true(all(expected %in% names(.get_result()$case_number_detail)))
})

test_that("case_number_detail satisfies the accounting invariant", {
  cnd <- .get_result()$case_number_detail
  expect_equal(
    cnd$n_original_df,
    cnd$n_filtered_out + length(cnd$error_arimax_skipped) + cnd$n_used_iarimax
  )
})

test_that("error_arimax_skipped is a character vector", {
  expect_type(.get_result()$case_number_detail$error_arimax_skipped, "character")
})

test_that("meta_analysis is a metafor rma object", {
  expect_s3_class(.get_result()$meta_analysis, "rma")
})

test_that("models is NULL when keep_models = FALSE", {
  expect_null(.get_result()$models)
})

test_that("models list has correct length and Arima class when keep_models = TRUE", {
  skip_on_cran()
  res      <- iarimax(dataframe = panel, y_series = "y", x_series = "x",
                      id_var = "id", timevar = "time", keep_models = TRUE)
  non_null <- Filter(Negate(is.null), res$models)
  # Total list length matches nrow(results_df) — includes NULL slots for failed fits
  expect_length(res$models, nrow(res$results_df))
  # Non-NULL count matches n_used_iarimax
  expect_equal(length(non_null), res$case_number_detail$n_used_iarimax)
  # All retained models carry the Arima class
  expect_true(all(sapply(non_null, function(m) "Arima" %in% class(m))))
})


# ── results_df columns ────────────────────────────────────────────────────────

test_that("results_df contains all required columns", {
  required <- c("id", "nAR", "nI", "nMA", "raw_cor",
                "n_valid", "n_params", "estimate_x", "std.error_x")
  expect_true(all(required %in% colnames(.get_result()$results_df)))
})

test_that("results_df id column is character type", {
  expect_type(.get_result()$results_df$id, "character")
})

test_that("results_df has one row per subject that passed filtering", {
  result <- .get_result()
  expect_equal(nrow(result$results_df), length(unique(panel$id)))
})

test_that("n_valid is non-negative for all fitted subjects", {
  valid_rows <- .get_result()$results_df
  valid_rows <- valid_rows[!is.na(valid_rows$n_valid), ]
  expect_true(all(valid_rows$n_valid >= 0))
})

test_that("verbose = TRUE produces same numeric results as verbose = FALSE", {
  skip_on_cran()
  result   <- .get_result()
  res_v    <- suppressMessages(
    iarimax(dataframe = panel, y_series = "y", x_series = "x",
            id_var = "id", timevar = "time", verbose = TRUE)
  )
  expect_equal(result$results_df$estimate_x,  res_v$results_df$estimate_x)
  expect_equal(result$results_df$std.error_x, res_v$results_df$std.error_x)
  expect_equal(result$results_df$raw_cor,     res_v$results_df$raw_cor)
})
