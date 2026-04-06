# Tests for i_screener()
# No model fitting — all tests run fast; no skip_on_cran() needed.

# ── Shared helper ─────────────────────────────────────────────────────────────
# Five subjects covering each quality failure mode:
#   "A" — passes all criteria on x and y
#   "B" — x has only 2 valid obs       (fails min_n when min_n >= 3)
#   "C" — x is nearly constant         (fails min_sd; also fails max_mode_pct)
#   "D" — x has 4/5 identical values   (fails max_mode_pct; passes min_sd)
#   "E" — y is nearly constant         (fails min_sd on y; x is fine)
#
# Known metric values (used in correctness tests):
#   A x: n=5, SD≈1.581, mode_pct=0.20
#   B x: n=2, SD≈1.414, mode_pct=0.50
#   C x: n=5, SD≈0.447, mode_pct=0.80   (mean=2.2: four 2s, one 3)
#   D x: n=5, SD≈0.894, mode_pct=0.80   (four 5s, one 3)
#   E y: n=5, SD≈0.447, mode_pct=0.80   (four 3s, one 4)

make_screen_df <- function() {
  data.frame(
    id = rep(c("A", "B", "C", "D", "E"), each = 5),
    x = c(
      1, 2, 4, 3, 5,       # A: varied
      1, 3, NA, NA, NA,    # B: only 2 valid
      2, 2, 2, 2, 3,       # C: nearly constant, SD ≈ 0.447
      5, 5, 5, 5, 3,       # D: mode_pct = 0.80, SD ≈ 0.894
      1, 2, 4, 3, 5        # E: varied (same as A)
    ),
    y = c(
      5, 3, 1, 4, 2,       # A: varied
      2, 4, 1, 5, 3,       # B: all 5 valid
      2, 4, 1, 5, 3,       # C: varied
      2, 4, 1, 5, 3,       # D: varied
      3, 3, 3, 3, 4        # E: nearly constant, SD ≈ 0.447
    ),
    stringsAsFactors = FALSE
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# Validation
# ══════════════════════════════════════════════════════════════════════════════

test_that("misspelled col triggers error naming the bad variable", {
  df <- make_screen_df()
  expect_error(
    i_screener(df, cols = "not_a_col", idvar = "id"),
    regexp = "not_a_col"
  )
})

test_that("misspelled idvar triggers error naming the bad variable", {
  df <- make_screen_df()
  expect_error(
    i_screener(df, cols = "x", idvar = "not_an_id"),
    regexp = "not_an_id"
  )
})

test_that("error message includes 'Cannot find required variables'", {
  df <- make_screen_df()
  expect_error(
    i_screener(df, cols = "oops", idvar = "id"),
    regexp = "Cannot find required variables"
  )
})

test_that("empty cols vector triggers error mentioning 'cols'", {
  df <- make_screen_df()
  expect_error(
    i_screener(df, cols = character(0), idvar = "id"),
    regexp = "cols"
  )
})

test_that("max_mode_pct = 0 triggers error", {
  df <- make_screen_df()
  expect_error(
    i_screener(df, cols = "x", idvar = "id", max_mode_pct = 0),
    regexp = "max_mode_pct"
  )
})

test_that("max_mode_pct > 1 triggers error", {
  df <- make_screen_df()
  expect_error(
    i_screener(df, cols = "x", idvar = "id", max_mode_pct = 1.1),
    regexp = "max_mode_pct"
  )
})

test_that("min_sd <= 0 triggers error", {
  df <- make_screen_df()
  expect_error(
    i_screener(df, cols = "x", idvar = "id", min_sd = 0),
    regexp = "min_sd"
  )
})

test_that("min_n_subject = NA_real_ triggers user-friendly error", {
  df <- make_screen_df()
  expect_error(
    i_screener(df, cols = "x", idvar = "id", min_n_subject = NA_real_),
    regexp = "min_n"
  )
})

test_that("min_n_subject = Inf triggers user-friendly error", {
  df <- make_screen_df()
  expect_error(
    i_screener(df, cols = "x", idvar = "id", min_n_subject = Inf),
    regexp = "min_n"
  )
})

test_that("min_sd = NA_real_ triggers user-friendly error", {
  df <- make_screen_df()
  expect_error(
    i_screener(df, cols = "x", idvar = "id", min_n_subject = 3, min_sd = NA_real_),
    regexp = "min_sd"
  )
})

test_that("idvar as a vector triggers error", {
  df <- make_screen_df()
  expect_error(
    i_screener(df, cols = "x", idvar = c("id", "id")),
    regexp = "idvar"
  )
})

test_that("min_sd = Inf triggers user-friendly error", {
  df <- make_screen_df()
  expect_error(
    i_screener(df, cols = "x", idvar = "id", min_sd = Inf),
    regexp = "min_sd"
  )
})

test_that("invalid mode triggers error", {
  df <- make_screen_df()
  expect_error(
    i_screener(df, cols = "x", idvar = "id", mode = "bad_mode"),
    regexp = "mode"
  )
})

test_that("invalid filter_type triggers error", {
  df <- make_screen_df()
  expect_error(
    i_screener(df, cols = "x", idvar = "id", filter_type = "bad_type"),
    regexp = "filter_type"
  )
})


# ══════════════════════════════════════════════════════════════════════════════
# Output structure
# ══════════════════════════════════════════════════════════════════════════════

test_that("mode=filter joint returns a data.frame", {
  df     <- make_screen_df()
  result <- i_screener(df, cols = "x", idvar = "id", min_n_subject = 3)
  expect_s3_class(result, "data.frame")
})

test_that("mode=filter joint preserves all original columns", {
  df     <- make_screen_df()
  result <- i_screener(df, cols = "x", idvar = "id", min_n_subject = 3)
  expect_true(all(c("id", "x", "y") %in% names(result)))
  expect_equal(ncol(result), ncol(df))
})

test_that("mode=filter joint removes rows of failing subjects", {
  df     <- make_screen_df()
  # B has n_valid = 2 for x; with min_n_subject = 3 it fails; others pass
  result <- i_screener(df, cols = "x", idvar = "id", min_n_subject = 3)
  expect_false("B" %in% result$id)
  expect_true(all(c("A", "C", "D", "E") %in% result$id))
  expect_equal(nrow(result), 20L)  # 4 subjects x 5 rows
})

test_that("mode=filter per_column does not remove rows", {
  df     <- make_screen_df()
  result <- i_screener(df, cols = "x", idvar = "id", min_n_subject = 3,
                     filter_type = "per_column")
  expect_equal(nrow(result), nrow(df))
  expect_equal(ncol(result), ncol(df))
})

test_that("mode=filter per_column sets failing col values to NA, other cols unchanged", {
  df <- make_screen_df()
  # B fails min_n on x; its y values should remain intact
  result <- i_screener(df, cols = "x", idvar = "id", min_n_subject = 3,
                     filter_type = "per_column")
  x_B <- result$x[result$id == "B"]
  y_B <- result$y[result$id == "B"]
  expect_true(all(is.na(x_B)))
  expect_false(any(is.na(y_B)))
})

test_that("mode=flag joint adds exactly one pass_overall column", {
  df     <- make_screen_df()
  result <- i_screener(df, cols = "x", idvar = "id", min_n_subject = 3, mode = "flag")
  expect_true("pass_overall" %in% names(result))
  expect_equal(ncol(result), ncol(df) + 1L)
  expect_equal(nrow(result), nrow(df))
  expect_type(result$pass_overall, "logical")
})

test_that("mode=flag per_column adds one <col>_pass column per variable", {
  df     <- make_screen_df()
  result <- i_screener(df, cols = c("x", "y"), idvar = "id", min_n_subject = 3,
                     mode = "flag", filter_type = "per_column")
  expect_true("x_pass" %in% names(result))
  expect_true("y_pass" %in% names(result))
  expect_false("pass_overall" %in% names(result))
  expect_equal(ncol(result), ncol(df) + 2L)
})

test_that("mode=report returns one row per subject", {
  df     <- make_screen_df()
  result <- i_screener(df, cols = c("x", "y"), idvar = "id", min_n_subject = 3,
                     mode = "report")
  expect_equal(nrow(result), 5L)  # 5 subjects
})

test_that("mode=report contains expected columns", {
  df     <- make_screen_df()
  result <- i_screener(df, cols = c("x", "y"), idvar = "id", min_n_subject = 3,
                     mode = "report")
  expected_cols <- c("id",
                     "x_n_valid", "y_n_valid",
                     "x_sd",      "y_sd",
                     "x_mode_pct","y_mode_pct",
                     "x_pass",   "y_pass",
                     "pass_overall")
  expect_true(all(expected_cols %in% names(result)))
})

test_that("output has no residual grouping", {
  df     <- make_screen_df()
  result <- i_screener(df, cols = "x", idvar = "id", min_n_subject = 3)
  expect_false(dplyr::is_grouped_df(result))
})

test_that("pre-existing grouping on input is handled without error", {
  df     <- dplyr::group_by(make_screen_df(), id)
  result <- i_screener(df, cols = "x", idvar = "id", min_n_subject = 3)
  expect_s3_class(result, "data.frame")
})

test_that("row order of passing subjects is preserved", {
  df     <- make_screen_df()
  result <- i_screener(df, cols = "x", idvar = "id", min_n_subject = 3)
  # All rows from B are removed; remaining rows keep original order
  df_expected <- df[df$id != "B", ]
  expect_equal(result$id,   df_expected$id)
  expect_equal(result$x,    df_expected$x,    ignore_attr = TRUE)
})


# ══════════════════════════════════════════════════════════════════════════════
# Statistical correctness
# ══════════════════════════════════════════════════════════════════════════════

test_that("min_n removes subjects below the threshold (joint)", {
  df     <- make_screen_df()
  # B: x has n_valid = 2; min_n_subject = 3 → B fails
  result <- i_screener(df, cols = "x", idvar = "id", min_n_subject = 3)
  expect_false("B" %in% result$id)
})

test_that("min_n keeps subjects at or above the threshold", {
  df     <- make_screen_df()
  result <- i_screener(df, cols = "x", idvar = "id", min_n_subject = 2)
  expect_true(all(c("A", "B", "C", "D", "E") %in% result$id))
})

test_that("min_sd removes low-variance subjects (joint)", {
  df     <- make_screen_df()
  # C: x SD ≈ 0.447 < 0.5 → C fails; D: SD ≈ 0.894 → passes
  result <- i_screener(df, cols = "x", idvar = "id",
                     min_n_subject = 2, min_sd = 0.5)
  expect_false("C" %in% result$id)
  expect_true("D"  %in% result$id)
})

test_that("min_sd removes subject with low SD on second col (joint)", {
  df     <- make_screen_df()
  # E: y SD ≈ 0.447 < 0.5 → E fails in joint mode with cols = c("x", "y")
  result <- i_screener(df, cols = c("x", "y"), idvar = "id",
                     min_n_subject = 2, min_sd = 0.5)
  expect_false("E" %in% result$id)
})

test_that("min_sd per_column sets NA only in failing column", {
  df <- make_screen_df()
  # C: x fails min_sd; y is fine → x values NA, y unchanged
  result <- i_screener(df, cols = c("x", "y"), idvar = "id",
                     min_n_subject = 2, min_sd = 0.5,
                     filter_type = "per_column")
  x_C <- result$x[result$id == "C"]
  y_C <- result$y[result$id == "C"]
  expect_true(all(is.na(x_C)))
  expect_false(any(is.na(y_C)))
})

test_that("max_mode_pct removes stuck responders (joint)", {
  df <- make_screen_df()
  # C: mode_pct = 0.80 > 0.75 → fails; D: mode_pct = 0.80 > 0.75 → fails
  result <- i_screener(df, cols = "x", idvar = "id",
                     min_n_subject = 2, max_mode_pct = 0.75)
  expect_false("C" %in% result$id)
  expect_false("D" %in% result$id)
  expect_true("A"  %in% result$id)
})

test_that("max_mode_pct at exactly the threshold value passes", {
  df <- make_screen_df()
  # C and D both have mode_pct = 0.80; max_mode_pct = 0.80 → they pass
  result <- i_screener(df, cols = "x", idvar = "id",
                     min_n_subject = 2, max_mode_pct = 0.80)
  expect_true("C" %in% result$id)
  expect_true("D" %in% result$id)
})

test_that("max_mode_pct per_column sets NA only in failing column", {
  df <- make_screen_df()
  # D: x fails max_mode_pct; D's y is fine
  result <- i_screener(df, cols = c("x", "y"), idvar = "id",
                     min_n_subject = 2, max_mode_pct = 0.75,
                     filter_type = "per_column")
  x_D <- result$x[result$id == "D"]
  y_D <- result$y[result$id == "D"]
  expect_true(all(is.na(x_D)))
  expect_false(any(is.na(y_D)))
})

test_that("joint mode removes subject failing on second col only", {
  df <- make_screen_df()
  # E: x passes min_sd, y fails min_sd → joint removes E; per_column keeps E rows
  result_joint <- i_screener(df, cols = c("x", "y"), idvar = "id",
                            min_n_subject = 2, min_sd = 0.5, filter_type = "joint")
  result_pc    <- i_screener(df, cols = c("x", "y"), idvar = "id",
                            min_n_subject = 2, min_sd = 0.5, filter_type = "per_column")
  expect_false("E" %in% result_joint$id)
  expect_true("E"  %in% result_pc$id)
})

test_that("flag joint: pass_overall is TRUE only for fully passing subjects", {
  df     <- make_screen_df()
  result <- i_screener(df, cols = c("x", "y"), idvar = "id",
                     min_n_subject = 2, min_sd = 0.5, mode = "flag")
  # A passes both; B, C, E fail; D passes x but D's y is fine
  expect_true( all(result$pass_overall[result$id == "A"]))
  expect_false(all(result$pass_overall[result$id == "C"]))
  expect_false(all(result$pass_overall[result$id == "E"]))
})

test_that("flag per_column: x_pass and y_pass reflect independent evaluation", {
  df     <- make_screen_df()
  result <- i_screener(df, cols = c("x", "y"), idvar = "id",
                     min_n_subject = 2, min_sd = 0.5,
                     mode = "flag", filter_type = "per_column")
  # E: x passes (SD ≈ 1.58 > 0.5), y fails (SD ≈ 0.447 < 0.5)
  x_pass_E <- unique(result$x_pass[result$id == "E"])
  y_pass_E <- unique(result$y_pass[result$id == "E"])
  expect_true(x_pass_E)
  expect_false(y_pass_E)
})

test_that("report returns correct n_valid values", {
  df     <- make_screen_df()
  result <- i_screener(df, cols = "x", idvar = "id", min_n_subject = 2, mode = "report")
  expect_equal(result$x_n_valid[result$id == "A"], 5L)
  expect_equal(result$x_n_valid[result$id == "B"], 2L)
})

test_that("report returns correct sd values", {
  df     <- make_screen_df()
  result <- i_screener(df, cols = "x", idvar = "id", min_n_subject = 2, mode = "report")
  expect_equal(result$x_sd[result$id == "A"],
               sd(c(1, 2, 4, 3, 5)), tolerance = 1e-10)
  expect_equal(result$x_sd[result$id == "C"],
               sd(c(2, 2, 2, 2, 3)), tolerance = 1e-10)
})

test_that("report returns correct mode_pct values", {
  df     <- make_screen_df()
  result <- i_screener(df, cols = "x", idvar = "id", min_n_subject = 2, mode = "report")
  # A: all unique → 1/5 = 0.2; C: four 2s → 4/5 = 0.8; D: four 5s → 0.8
  expect_equal(result$x_mode_pct[result$id == "A"], 0.2, tolerance = 1e-10)
  expect_equal(result$x_mode_pct[result$id == "C"], 0.8, tolerance = 1e-10)
  expect_equal(result$x_mode_pct[result$id == "D"], 0.8, tolerance = 1e-10)
})

test_that("report pass_overall is FALSE when any col fails", {
  df     <- make_screen_df()
  result <- i_screener(df, cols = c("x", "y"), idvar = "id",
                     min_n_subject = 2, min_sd = 0.5, mode = "report")
  # E passes x but fails y → pass_overall must be FALSE
  expect_false(result$pass_overall[result$id == "E"])
})

test_that("no criteria beyond min_n returns all subjects with n_valid >= min_n", {
  df     <- make_screen_df()
  result <- i_screener(df, cols = "x", idvar = "id", min_n_subject = 3)
  # Only B (n=2) is removed; C, D, E all have n=5
  expect_true(all(c("A", "C", "D", "E") %in% result$id))
  expect_false("B" %in% result$id)
})


# ══════════════════════════════════════════════════════════════════════════════
# Verbose
# ══════════════════════════════════════════════════════════════════════════════

test_that("verbose = TRUE emits messages", {
  df <- make_screen_df()
  suppressMessages(
    expect_message(
      i_screener(df, cols = "x", idvar = "id", min_n_subject = 3, verbose = TRUE)
    )
  )
})

test_that("verbose = FALSE emits no messages", {
  df <- make_screen_df()
  expect_no_message(
    i_screener(df, cols = "x", idvar = "id", min_n_subject = 3, verbose = FALSE)
  )
})

test_that("verbose message mentions 'i_screener' and subject counts", {
  df <- make_screen_df()
  msgs <- character(0)
  withCallingHandlers(
    i_screener(df, cols = "x", idvar = "id", min_n_subject = 3, verbose = TRUE),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  combined <- paste(msgs, collapse = " ")
  expect_match(combined, "i_screener",  ignore.case = TRUE)
  expect_match(combined, "min_n",       ignore.case = TRUE)
  expect_match(combined, "[0-9]")  # at least one count in the output
})


# ══════════════════════════════════════════════════════════════════════════════
# Column name collision guard
# ══════════════════════════════════════════════════════════════════════════════

test_that("flag joint errors if df already has pass_overall column", {
  df              <- make_screen_df()
  df$pass_overall  <- TRUE
  expect_error(
    i_screener(df, cols = "x", idvar = "id", min_n_subject = 3, mode = "flag"),
    regexp = "pass_overall"
  )
})

test_that("flag per_column errors if df already has a <col>_pass column", {
  df       <- make_screen_df()
  df$x_pass <- TRUE
  expect_error(
    i_screener(df, cols = "x", idvar = "id", min_n_subject = 3,
             mode = "flag", filter_type = "per_column"),
    regexp = "x_pass"
  )
})

test_that("filter per_column errors if df already has a <col>_pass column", {
  df        <- make_screen_df()
  df$x_pass <- TRUE
  expect_error(
    i_screener(df, cols = "x", idvar = "id", min_n_subject = 3,
             mode = "filter", filter_type = "per_column"),
    regexp = "x_pass"
  )
})


# ══════════════════════════════════════════════════════════════════════════════
# Additional correctness edge cases
# ══════════════════════════════════════════════════════════════════════════════

test_that("mode=filter returns all rows when all subjects pass", {
  df     <- make_screen_df()
  result <- i_screener(df, cols = "x", idvar = "id", min_n_subject = 1)
  expect_equal(nrow(result), nrow(df))
  expect_equal(result$id, df$id)
})

test_that("max_mode_pct = 1.0 does not remove subjects with mode_pct < 1.0", {
  df     <- make_screen_df()
  # C has mode_pct = 0.80 ≤ 1.0 → passes; everyone with n ≥ 2 should pass
  result <- i_screener(df, cols = "x", idvar = "id",
                     min_n_subject = 2, max_mode_pct = 1.0)
  expect_true("C" %in% result$id)
  expect_true("D" %in% result$id)
})

test_that("subject with n_valid = 1 fails min_sd due to NA sd", {
  df <- data.frame(
    id = c(rep("A", 5), rep("B", 5)),
    x  = c(1, 2, 3, 4, 5,          # A: n=5, SD > 0
            5, NA, NA, NA, NA),     # B: n=1, SD = NA → fails min_sd
    stringsAsFactors = FALSE
  )
  result <- i_screener(df, cols = "x", idvar = "id",
                     min_n_subject = 1, min_sd = 0.5, mode = "report")
  expect_false(result$x_pass[result$id == "B"])
  expect_true(result$x_pass[result$id == "A"])
})

test_that("per_column filter handles different subjects failing different cols independently", {
  df <- data.frame(
    id = rep(c("A", "B"), each = 5),
    x  = c(1, 2, 3, 4, 5,    # A: SD > 0.5 → passes
            3, 3, 3, 3, 3),   # B: constant → SD = 0 → fails min_sd on x
    y  = c(4, 4, 4, 4, 4,    # A: constant → SD = 0 → fails min_sd on y
            1, 2, 3, 4, 5),   # B: varied → passes
    stringsAsFactors = FALSE
  )
  result <- i_screener(df, cols = c("x", "y"), idvar = "id",
                     min_n_subject = 2, min_sd = 0.5,
                     filter_type = "per_column")
  # A: y is NA; x is intact
  expect_true(all(is.na(result$y[result$id == "A"])))
  expect_false(any(is.na(result$x[result$id == "A"])))
  # B: x is NA; y is intact
  expect_true(all(is.na(result$x[result$id == "B"])))
  expect_false(any(is.na(result$y[result$id == "B"])))
})

test_that("report column order groups by metric type across all cols", {
  df     <- make_screen_df()
  result <- i_screener(df, cols = c("x", "y"), idvar = "id",
                     min_n_subject = 2, mode = "report")
  expected_order <- c("id",
                      "x_n_valid", "y_n_valid",
                      "x_sd",      "y_sd",
                      "x_mode_pct","y_mode_pct",
                      "x_pass",    "y_pass",
                      "pass_overall")
  expect_equal(names(result), expected_order)
})
