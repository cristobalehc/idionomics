# Tests for i_detrender()
# No skip_on_cran() needed — only lm() fitting, which is fast.

# ── Shared helper ──────────────────────────────────────────────────────────────
# Five subjects covering the main filter cases:
#   "valid"      — 25 obs, normal variance, detrends cleanly
#   "too_few"    — 5 obs, below default min_n_subject = 20
#   "no_var"     — 25 obs, constant x (pre-detrend var = 0)
#   "pure_trend" — 25 obs, perfect linear x (post-detrend var = 0)
#   "all_na"     — 25 obs, all-NA x

make_det_df <- function() {
  set.seed(1)
  rbind(
    data.frame(id = "valid",      time = seq_len(25), x = rnorm(25),               stringsAsFactors = FALSE),
    data.frame(id = "too_few",    time = seq_len(5),  x = rnorm(5),                stringsAsFactors = FALSE),
    data.frame(id = "no_var",     time = seq_len(25), x = rep(3.0, 25),            stringsAsFactors = FALSE),
    data.frame(id = "pure_trend", time = seq_len(25), x = as.numeric(seq_len(25)), stringsAsFactors = FALSE),
    data.frame(id = "all_na",     time = seq_len(25), x = rep(NA_real_, 25),       stringsAsFactors = FALSE)
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# Validation
# ══════════════════════════════════════════════════════════════════════════════

test_that("misspelled col triggers error naming the bad variable", {
  df <- make_det_df()
  expect_error(
    i_detrender(df, cols = "not_a_col", id_var = "id", timevar = "time"),
    regexp = "not_a_col"
  )
})

test_that("misspelled id_var triggers error naming the bad variable", {
  df <- make_det_df()
  expect_error(
    i_detrender(df, cols = "x", id_var = "not_an_id", timevar = "time"),
    regexp = "not_an_id"
  )
})

test_that("misspelled timevar triggers error naming the bad variable", {
  df <- make_det_df()
  expect_error(
    i_detrender(df, cols = "x", id_var = "id", timevar = "not_a_time"),
    regexp = "not_a_time"
  )
})

test_that("error message includes 'Cannot find required variables'", {
  df <- make_det_df()
  expect_error(
    i_detrender(df, cols = "oops", id_var = "id", timevar = "time"),
    regexp = "Cannot find required variables"
  )
})

test_that("empty cols vector triggers error", {
  df <- make_det_df()
  expect_error(
    i_detrender(df, cols = character(0), id_var = "id", timevar = "time"),
    regexp = "cols"
  )
})

test_that("character column in cols triggers error naming the column", {
  df <- make_det_df()
  df$x <- as.character(df$x)
  expect_error(
    i_detrender(df, cols = "x", id_var = "id", timevar = "time"),
    regexp = "must be numeric"
  )
})

test_that("non-numeric cols error reports the actual class", {
  df <- make_det_df()
  df$x <- factor(df$x)
  expect_error(
    i_detrender(df, cols = "x", id_var = "id", timevar = "time"),
    regexp = "factor"
  )
})

test_that("id_var as a vector triggers error", {
  df <- make_det_df()
  expect_error(
    i_detrender(df, cols = "x", id_var = c("id", "id"), timevar = "time"),
    regexp = "id_var"
  )
})

test_that("timevar as a vector triggers error", {
  df <- make_det_df()
  expect_error(
    i_detrender(df, cols = "x", id_var = "id", timevar = c("time", "time")),
    regexp = "timevar"
  )
})

test_that("min_n_subject = 0 triggers upfront error", {
  df <- make_det_df()
  expect_error(
    i_detrender(df, cols = "x", id_var = "id", timevar = "time", min_n_subject = 0),
    regexp = "min_n_subject"
  )
})

test_that("minvar = -1 triggers upfront error", {
  df <- make_det_df()
  expect_error(
    i_detrender(df, cols = "x", id_var = "id", timevar = "time", minvar = -1),
    regexp = "minvar"
  )
})

test_that("non-numeric timevar triggers error naming the class", {
  df <- make_det_df()
  df$time <- as.character(df$time)
  expect_error(
    i_detrender(df, cols = "x", id_var = "id", timevar = "time"),
    regexp = "must be numeric"
  )
})

test_that("Date timevar triggers error suggesting numeric conversion", {
  df <- make_det_df()
  df$time <- as.Date("2024-01-01") + df$time
  expect_error(
    i_detrender(df, cols = "x", id_var = "id", timevar = "time"),
    regexp = "must be numeric"
  )
})

test_that("single NA in timevar triggers error reporting the count", {
  df         <- make_det_df()
  df$time[1] <- NA
  expect_error(
    i_detrender(df, cols = "x", id_var = "id", timevar = "time"),
    regexp = "1 row"
  )
})

test_that("NA in timevar error names the offending variable", {
  df         <- make_det_df()
  df$time[1] <- NA
  expect_error(
    i_detrender(df, cols = "x", id_var = "id", timevar = "time"),
    regexp = "'time'"
  )
})


# ══════════════════════════════════════════════════════════════════════════════
# Output structure
# ══════════════════════════════════════════════════════════════════════════════

test_that("returns a data.frame", {
  df     <- make_det_df()
  result <- i_detrender(df, cols = "x", id_var = "id", timevar = "time")
  expect_s3_class(result, "data.frame")
})

test_that("append = TRUE preserves all original columns", {
  df     <- make_det_df()
  result <- i_detrender(df, cols = "x", id_var = "id", timevar = "time", append = TRUE)
  expect_true(all(c("id", "time", "x") %in% names(result)))
})

test_that("append = TRUE adds exactly one _dt column for one input col", {
  df     <- make_det_df()
  result <- i_detrender(df, cols = "x", id_var = "id", timevar = "time", append = TRUE)
  expect_true("x_dt" %in% names(result))
  expect_equal(ncol(result), ncol(df) + 1L)
})

test_that("append = TRUE preserves row count", {
  df     <- make_det_df()
  result <- i_detrender(df, cols = "x", id_var = "id", timevar = "time", append = TRUE)
  expect_equal(nrow(result), nrow(df))
})

test_that("append = FALSE returns only id_var, timevar, and _dt columns", {
  df     <- make_det_df()
  result <- i_detrender(df, cols = "x", id_var = "id", timevar = "time", append = FALSE)
  expect_equal(sort(names(result)), sort(c("id", "time", "x_dt")))
})

test_that("append = FALSE preserves row count", {
  df     <- make_det_df()
  result <- i_detrender(df, cols = "x", id_var = "id", timevar = "time", append = FALSE)
  expect_equal(nrow(result), nrow(df))
})

test_that("multiple cols produce a _dt column for each", {
  set.seed(2)
  df     <- make_det_df()
  df$x2  <- c(rnorm(25), rnorm(5), rep(3.0, 25), as.numeric(seq_len(25)), rep(NA_real_, 25))
  result <- i_detrender(df, cols = c("x", "x2"), id_var = "id", timevar = "time")
  expect_true("x_dt"  %in% names(result))
  expect_true("x2_dt" %in% names(result))
  expect_equal(ncol(result), ncol(df) + 2L)
})

test_that("row order of input is preserved in output", {
  df     <- make_det_df()
  result <- i_detrender(df, cols = "x", id_var = "id", timevar = "time")
  expect_equal(result$id, df$id)
})

test_that("pre-existing grouping is removed without error", {
  df     <- dplyr::group_by(make_det_df(), id)
  result <- i_detrender(df, cols = "x", id_var = "id", timevar = "time")
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), nrow(make_det_df()))
})

test_that("output df has no residual grouping", {
  df     <- make_det_df()
  result <- i_detrender(df, cols = "x", id_var = "id", timevar = "time")
  expect_false(dplyr::is_grouped_df(result))
})


# ══════════════════════════════════════════════════════════════════════════════
# Per-subject filtering
# ══════════════════════════════════════════════════════════════════════════════

test_that("all-NA column within a subject produces all NA in _dt", {
  df     <- make_det_df()
  result <- i_detrender(df, cols = "x", id_var = "id", timevar = "time")
  expect_true(all(is.na(result$x_dt[result$id == "all_na"])))
})

test_that("subject below min_n_subject produces all NA in _dt", {
  df     <- make_det_df()   # "too_few" has 5 obs, default min_n_subject = 20
  result <- i_detrender(df, cols = "x", id_var = "id", timevar = "time")
  expect_true(all(is.na(result$x_dt[result$id == "too_few"])))
})

test_that("subject with pre-detrend variance < minvar produces all NA in _dt", {
  df     <- make_det_df()   # "no_var" has constant x, var = 0 < 0.01
  result <- i_detrender(df, cols = "x", id_var = "id", timevar = "time")
  expect_true(all(is.na(result$x_dt[result$id == "no_var"])))
})

test_that("subject with post-detrend variance < minvar produces all NA in _dt", {
  df     <- make_det_df()   # "pure_trend": seq(1,25) ~ seq(1,25) → residuals = 0
  result <- i_detrender(df, cols = "x", id_var = "id", timevar = "time")
  expect_true(all(is.na(result$x_dt[result$id == "pure_trend"])))
})

test_that("valid subject passes all filters and produces non-NA _dt values", {
  df     <- make_det_df()
  result <- i_detrender(df, cols = "x", id_var = "id", timevar = "time")
  expect_false(any(is.na(result$x_dt[result$id == "valid"])))
})

test_that("filtering is per-column: one col can be NA while another is valid for the same subject", {
  set.seed(3)
  # Subject with good x1 (normal variance) but constant x2 (zero variance)
  df <- data.frame(
    id   = rep("A", 25),
    time = seq_len(25),
    x1   = rnorm(25),
    x2   = rep(5.0, 25),
    stringsAsFactors = FALSE
  )
  result <- i_detrender(df, cols = c("x1", "x2"), id_var = "id", timevar = "time")
  expect_false(any(is.na(result$x1_dt)))   # x1 detrends cleanly
  expect_true(all(is.na(result$x2_dt)))    # x2 is constant → NA
})

test_that("minvar non-default threshold is respected for pre-detrend variance", {
  set.seed(10)
  df <- rbind(
    data.frame(id = "low_var", time = seq_len(25), x = rnorm(25, sd = 0.08), stringsAsFactors = FALSE),
    data.frame(id = "ok",      time = seq_len(25), x = rnorm(25),            stringsAsFactors = FALSE)
  )
  # sd ≈ 0.08 → var ≈ 0.0064 < default minvar 0.01 → should fail
  result_strict <- i_detrender(df, cols = "x", id_var = "id", timevar = "time", minvar = 0.01)
  expect_true(all(is.na(result_strict$x_dt[result_strict$id == "low_var"])))

  # With minvar = 0.001: same subject passes
  result_loose <- i_detrender(df, cols = "x", id_var = "id", timevar = "time", minvar = 0.001)
  expect_false(any(is.na(result_loose$x_dt[result_loose$id == "low_var"])))
})

test_that("append = FALSE with multiple cols returns id_var, timevar, and all _dt columns", {
  set.seed(7)
  df    <- make_det_df()
  df$x2 <- c(rnorm(25), rnorm(5), rep(3.0, 25), as.numeric(seq_len(25)), rep(NA_real_, 25))
  result <- i_detrender(df, cols = c("x", "x2"), id_var = "id", timevar = "time", append = FALSE)
  expect_equal(sort(names(result)), sort(c("id", "time", "x_dt", "x2_dt")))
})

test_that("min_n_subject threshold is respected with a non-default value", {
  set.seed(4)
  df <- rbind(
    data.frame(id = "borderline", time = seq_len(10), x = rnorm(10), stringsAsFactors = FALSE),
    data.frame(id = "ok",         time = seq_len(25), x = rnorm(25), stringsAsFactors = FALSE)
  )
  # With min_n_subject = 15: "borderline" (10 obs) should fail
  result_strict <- i_detrender(df, cols = "x", id_var = "id", timevar = "time", min_n_subject = 15)
  expect_true(all(is.na(result_strict$x_dt[result_strict$id == "borderline"])))

  # With min_n_subject = 5: "borderline" (10 obs) should pass
  result_loose <- i_detrender(df, cols = "x", id_var = "id", timevar = "time", min_n_subject = 5)
  expect_false(any(is.na(result_loose$x_dt[result_loose$id == "borderline"])))
})


# ══════════════════════════════════════════════════════════════════════════════
# Statistical correctness
# ══════════════════════════════════════════════════════════════════════════════

test_that("residuals of valid subject have mean close to 0", {
  df     <- make_det_df()
  result <- i_detrender(df, cols = "x", id_var = "id", timevar = "time")
  x_dt_valid <- result$x_dt[result$id == "valid"]
  expect_equal(mean(x_dt_valid, na.rm = TRUE), 0, tolerance = 1e-10)
})

test_that("residuals match manual lm() calculation for valid subject", {
  df     <- make_det_df()
  valid  <- df[df$id == "valid", ]
  manual <- stats::residuals(stats::lm(x ~ time, data = valid, na.action = stats::na.exclude))

  result <- i_detrender(df, cols = "x", id_var = "id", timevar = "time")
  expect_equal(
    as.numeric(result$x_dt[result$id == "valid"]),
    as.numeric(manual),
    tolerance = 1e-10
  )
})

test_that("detrending one subject does not alter another subject's residuals", {
  set.seed(5)
  df <- rbind(
    data.frame(id = "A", time = seq_len(25), x = rnorm(25), stringsAsFactors = FALSE),
    data.frame(id = "B", time = seq_len(25), x = rnorm(25), stringsAsFactors = FALSE)
  )
  result   <- i_detrender(df, cols = "x", id_var = "id", timevar = "time")
  manual_B <- stats::residuals(
    stats::lm(x ~ time, data = df[df$id == "B", ], na.action = stats::na.exclude)
  )
  expect_equal(
    as.numeric(result$x_dt[result$id == "B"]),
    as.numeric(manual_B),
    tolerance = 1e-10
  )
})

test_that("NA values within an otherwise-varying series are preserved as NA in _dt", {
  set.seed(6)
  df <- data.frame(
    id   = rep("A", 25),
    time = seq_len(25),
    x    = c(NA, rnorm(24)),   # first obs is NA, rest are normal
    stringsAsFactors = FALSE
  )
  result <- i_detrender(df, cols = "x", id_var = "id", timevar = "time")
  expect_true(is.na(result$x_dt[1]))          # NA position preserved
  expect_false(any(is.na(result$x_dt[-1])))   # other positions filled
})

test_that("scattered interior NAs are preserved at correct positions in _dt", {
  set.seed(99)
  x_raw        <- rnorm(25)
  x_raw[c(5, 12, 20)] <- NA
  df <- data.frame(id = "A", time = seq_len(25), x = x_raw,
                   stringsAsFactors = FALSE)
  result <- i_detrender(df, cols = "x", id_var = "id", timevar = "time")
  expect_true(all(is.na(result$x_dt[c(5, 12, 20)])))
  expect_false(any(is.na(result$x_dt[-c(5, 12, 20)])))
})

test_that("all-NA subject produces all-NA output regardless of min_n_subject threshold", {
  df <- data.frame(id = "allna", time = seq_len(25),
                   x  = rep(NA_real_, 25), stringsAsFactors = FALSE)
  result <- i_detrender(df, cols = "x", id_var = "id", timevar = "time",
                        min_n_subject = 1)
  expect_true(all(is.na(result$x_dt)))
})

test_that("append = FALSE preserves rows for filtered-out subjects as NA", {
  df     <- make_det_df()
  result <- i_detrender(df, cols = "x", id_var = "id", timevar = "time",
                        append = FALSE)
  expect_equal(nrow(result), nrow(df))
  expect_true(all(is.na(result$x_dt[result$id == "no_var"])))
})


# ══════════════════════════════════════════════════════════════════════════════
# Verbose
# ══════════════════════════════════════════════════════════════════════════════

test_that("verbose = TRUE emits messages", {
  df   <- make_det_df()
  msgs <- capture_messages(
    i_detrender(df, cols = "x", id_var = "id", timevar = "time", verbose = TRUE)
  )
  expect_gt(length(msgs), 0)
})

test_that("verbose = FALSE emits no messages", {
  df <- make_det_df()
  expect_no_message(
    i_detrender(df, cols = "x", id_var = "id", timevar = "time", verbose = FALSE)
  )
})

test_that("verbose message mentions detrend", {
  df   <- make_det_df()
  msgs <- capture_messages(
    i_detrender(df, cols = "x", id_var = "id", timevar = "time", verbose = TRUE)
  )
  expect_true(any(grepl("(?i)detrend", msgs, perl = TRUE)))
})
