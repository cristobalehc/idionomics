# Tests for pmstandardize()
# No model fitting — all tests run fast; no skip_on_cran() needed.

# ── Shared helper ─────────────────────────────────────────────────────────────
# Three subjects covering the main edge cases:
#   "A" — normal variance (x and y)
#   "B" — constant x (zero variance), normal y
#   "C" — all-NA x, normal y

make_pms_df <- function() {
  data.frame(
    id = rep(c("A", "B", "C"), each = 5),
    x  = c( 1,  2,  3,  4,  5,      # A: mean = 3, sd = sqrt(2.5)
            10, 10, 10, 10, 10,      # B: constant
            NA, NA, NA, NA, NA),     # C: all NA
    y  = c( 5,  4,  3,  2,  1,      # A: normal
             1,  2,  3,  4,  5,      # B: normal
             2,  4,  6,  8, 10),     # C: normal
    stringsAsFactors = FALSE
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# Validation
# ══════════════════════════════════════════════════════════════════════════════

test_that("misspelled col triggers error naming the bad variable", {
  df <- make_pms_df()
  expect_error(
    pmstandardize(df, cols = "not_a_col", idvar = "id"),
    regexp = "not_a_col"
  )
})

test_that("misspelled idvar triggers error naming the bad variable", {
  df <- make_pms_df()
  expect_error(
    pmstandardize(df, cols = "x", idvar = "not_an_id"),
    regexp = "not_an_id"
  )
})

test_that("error message includes 'Cannot find required variables'", {
  df <- make_pms_df()
  expect_error(
    pmstandardize(df, cols = "oops", idvar = "id"),
    regexp = "Cannot find required variables"
  )
})

test_that("empty cols vector triggers error", {
  df <- make_pms_df()
  expect_error(
    pmstandardize(df, cols = character(0), idvar = "id"),
    regexp = "cols"
  )
})


# ══════════════════════════════════════════════════════════════════════════════
# Output structure
# ══════════════════════════════════════════════════════════════════════════════

test_that("returns a data.frame", {
  df     <- make_pms_df()
  result <- pmstandardize(df, cols = "x", idvar = "id")
  expect_s3_class(result, "data.frame")
})

test_that("append = TRUE preserves all original columns", {
  df     <- make_pms_df()
  result <- pmstandardize(df, cols = "x", idvar = "id", append = TRUE)
  expect_true(all(c("id", "x", "y") %in% names(result)))
})

test_that("append = TRUE adds exactly one _psd column for one input col", {
  df     <- make_pms_df()
  result <- pmstandardize(df, cols = "x", idvar = "id", append = TRUE)
  expect_true("x_psd" %in% names(result))
  # Total columns = original 3 + 1 new
  expect_equal(ncol(result), ncol(df) + 1L)
})

test_that("append = TRUE preserves row count", {
  df     <- make_pms_df()
  result <- pmstandardize(df, cols = "x", idvar = "id", append = TRUE)
  expect_equal(nrow(result), nrow(df))
})

test_that("append = FALSE returns only idvar and _psd columns", {
  df     <- make_pms_df()
  result <- pmstandardize(df, cols = "x", idvar = "id", append = FALSE)
  expect_equal(sort(names(result)), sort(c("id", "x_psd")))
})

test_that("append = FALSE preserves row count", {
  df     <- make_pms_df()
  result <- pmstandardize(df, cols = "x", idvar = "id", append = FALSE)
  expect_equal(nrow(result), nrow(df))
})

test_that("multiple cols produce a _psd column for each", {
  df     <- make_pms_df()
  result <- pmstandardize(df, cols = c("x", "y"), idvar = "id")
  expect_true("x_psd" %in% names(result))
  expect_true("y_psd" %in% names(result))
  expect_equal(ncol(result), ncol(df) + 2L)
})

test_that("row order of input is preserved in output", {
  df     <- make_pms_df()
  result <- pmstandardize(df, cols = "x", idvar = "id")
  expect_equal(result$id, df$id)
})

test_that("pre-existing grouping on input df is handled without error", {
  df     <- dplyr::group_by(make_pms_df(), id)
  result <- pmstandardize(df, cols = "x", idvar = "id")
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 15L)
})

test_that("output df has no residual grouping", {
  df     <- make_pms_df()
  result <- pmstandardize(df, cols = "x", idvar = "id")
  expect_false(dplyr::is_grouped_df(result))
})


# ══════════════════════════════════════════════════════════════════════════════
# Statistical correctness
# ══════════════════════════════════════════════════════════════════════════════

test_that("within-person mean of z-scores is 0 for normal-variance subjects", {
  df     <- make_pms_df()
  result <- pmstandardize(df, cols = "x", idvar = "id")

  # Subject A has normal variance — within-person mean must be 0
  x_psd_A <- result$x_psd[result$id == "A"]
  expect_equal(mean(x_psd_A, na.rm = TRUE), 0, tolerance = 1e-10)
})

test_that("within-person SD of z-scores is 1 for normal-variance subjects", {
  df     <- make_pms_df()
  result <- pmstandardize(df, cols = "x", idvar = "id")

  x_psd_A <- result$x_psd[result$id == "A"]
  expect_equal(sd(x_psd_A, na.rm = TRUE), 1, tolerance = 1e-10)
})

test_that("z-score values match manual calculation", {
  df     <- make_pms_df()
  result <- pmstandardize(df, cols = "x", idvar = "id")

  x_A        <- df$x[df$id == "A"]
  expected_z <- (x_A - mean(x_A)) / sd(x_A)
  x_psd_A    <- result$x_psd[result$id == "A"]

  expect_equal(x_psd_A, expected_z, tolerance = 1e-10)
})

test_that("constant within-person series returns 0 at observed positions, NA elsewhere", {
  df <- data.frame(
    id = rep(c("A", "B"), each = 5),
    x  = c(1, 2, 3, 4, 5,       # A: normal
           10, 10, NA, 10, 10),  # B: constant with one NA
    stringsAsFactors = FALSE
  )
  result  <- pmstandardize(df, cols = "x", idvar = "id")
  x_psd_B <- result$x_psd[result$id == "B"]

  expect_true(all(x_psd_B[c(1, 2, 4, 5)] == 0))   # observed constant values → 0
  expect_true(is.na(x_psd_B[3]))                    # NA position preserved
})

test_that("all-NA within-person series returns NA", {
  df     <- make_pms_df()
  result <- pmstandardize(df, cols = "x", idvar = "id")

  x_psd_C <- result$x_psd[result$id == "C"]
  expect_true(all(is.na(x_psd_C)))
})

test_that("person with a single non-NA value returns 0 at non-NA position, NA elsewhere", {
  df <- data.frame(
    id = rep(c("A", "B"), each = 5),
    x  = c(1, 2, 3, 4, 5,      # A: normal
           NA, NA, 7, NA, NA),  # B: exactly one non-NA
    stringsAsFactors = FALSE
  )
  result  <- pmstandardize(df, cols = "x", idvar = "id")
  x_psd_B <- result$x_psd[result$id == "B"]

  expect_equal(x_psd_B[3], 0)                          # the one observed value → 0
  expect_true(all(is.na(x_psd_B[c(1, 2, 4, 5)])))     # NA positions preserved
})

test_that("NA values within an otherwise-varying series are preserved as NA in output", {
  df <- data.frame(
    id = rep("A", 6),
    x  = c(1, NA, 3, 4, 5, 6),
    stringsAsFactors = FALSE
  )
  result <- pmstandardize(df, cols = "x", idvar = "id")

  expect_true(is.na(result$x_psd[2]))           # NA position preserved
  expect_false(any(is.na(result$x_psd[-2])))    # other positions filled
})

test_that("NA values do not contaminate within-person mean and SD calculation", {
  df <- data.frame(
    id = rep("A", 6),
    x  = c(1, NA, 3, 4, 5, 6),
    stringsAsFactors = FALSE
  )
  result <- pmstandardize(df, cols = "x", idvar = "id")

  non_na_x   <- df$x[!is.na(df$x)]
  expected_z <- (non_na_x - mean(non_na_x)) / sd(non_na_x)

  expect_equal(result$x_psd[!is.na(result$x_psd)], expected_z, tolerance = 1e-10)
})

test_that("each subject is standardized independently (between-person means differ)", {
  # Two subjects with very different scales — output should both have mean 0
  df <- data.frame(
    id = rep(c("A", "B"), each = 5),
    x  = c(1, 2, 3, 4, 5,          # A: scale ~1
           100, 200, 300, 400, 500), # B: scale ~100
    stringsAsFactors = FALSE
  )
  result <- pmstandardize(df, cols = "x", idvar = "id")

  mean_A <- mean(result$x_psd[result$id == "A"])
  mean_B <- mean(result$x_psd[result$id == "B"])

  expect_equal(mean_A, 0, tolerance = 1e-10)
  expect_equal(mean_B, 0, tolerance = 1e-10)

  sd_A <- sd(result$x_psd[result$id == "A"])
  sd_B <- sd(result$x_psd[result$id == "B"])

  expect_equal(sd_A, 1, tolerance = 1e-10)
  expect_equal(sd_B, 1, tolerance = 1e-10)
})

test_that("standardizing y does not alter x values and vice versa", {
  df     <- make_pms_df()
  result <- pmstandardize(df, cols = c("x", "y"), idvar = "id")

  expect_equal(result$x, df$x)
  expect_equal(result$y, df$y)
})


# ══════════════════════════════════════════════════════════════════════════════
# Verbose
# ══════════════════════════════════════════════════════════════════════════════

test_that("verbose = TRUE emits messages", {
  df <- make_pms_df()
  suppressMessages(
    expect_message(
      pmstandardize(df, cols = "x", idvar = "id", verbose = TRUE)
    )
  )
})

test_that("verbose = FALSE emits no messages", {
  df <- make_pms_df()
  expect_no_message(
    pmstandardize(df, cols = "x", idvar = "id", verbose = FALSE)
  )
})

test_that("verbose message mentions z-score or standardiz", {
  df <- make_pms_df()
  suppressMessages(
    expect_message(
      pmstandardize(df, cols = "x", idvar = "id", verbose = TRUE),
      regexp = "(?i)z.score|standardiz"
    )
  )
})
