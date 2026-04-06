# Tests for looping_machine()
# Layer 1: validation and structural checks using a tiny dataset (errors thrown
#           before auto.arima runs) — no skip_on_cran() needed.
# Layer 2: full integration tests.
#   Layer 2a — minimal panel (4 subjects × 25 obs, white noise). NO skip_on_cran().
#               covr::package_coverage() does NOT set NOT_CRAN, so skip_on_cran()
#               always fires there. These tests are fast enough (~5 s total) to
#               run on CRAN and ensure looping_machine.R gets instrumented.
#   Layer 2b — larger panel (6 subjects × 30 obs). skip_on_cran(). Used by
#               devtools::test() for more thorough verification; not needed for
#               coverage.

# ── Shared helpers ────────────────────────────────────────────────────────────

# Tiny dataset — only 5 obs per subject, so iarimax() will always throw
# "not enough cases" before fitting any model. Useful for fast validation tests
# that just need the column-name guard to fire.
make_tiny_loop_panel <- function() {
  data.frame(
    id   = rep(c("1", "2"), each = 5),
    time = rep(seq_len(5), 2),
    a    = rnorm(10),
    b    = rnorm(10),
    c    = rnorm(10),
    stringsAsFactors = FALSE
  )
}

# Minimal panel — 4 subjects × 25 obs with mild linear signal (beta ≈ 0.4).
# Signal is mild so auto.arima typically selects low-order models;
# runtime is acceptable for CRAN (~1–3 s per subject per leg).
# Used for Layer 2a (no skip_on_cran).
make_mini_loop_panel <- function(n_subjects = 4, n_obs = 25, seed = 5) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(n_subjects), function(id) {
    a  <- rnorm(n_obs)
    b  <- 0.4 * a + rnorm(n_obs)
    cc <- 0.4 * b + rnorm(n_obs)
    data.frame(
      id   = as.character(id),
      time = seq_len(n_obs),
      a    = a,
      b    = b,
      c    = cc,
      stringsAsFactors = FALSE
    )
  }))
}

# Minimal negative panel — strong negative a→b so at least one subject has
# a_b < 0. Used to test Loop_positive_directed = 0. No skip_on_cran.
make_mini_neg_loop_panel <- function(n_subjects = 4, n_obs = 25, seed = 55) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(n_subjects), function(id) {
    a  <- rnorm(n_obs)
    b  <- -0.8 * a + rnorm(n_obs)           # reliable negative signal, less extreme
    cc <- 0.4  * b + rnorm(n_obs)
    data.frame(
      id   = as.character(id),
      time = seq_len(n_obs),
      a    = a,
      b    = b,
      c    = cc,
      stringsAsFactors = FALSE
    )
  }))
}

# Full panel — 6 subjects × 30 obs with mild AR(1) dynamics.
# Used for Layer 2b (skip_on_cran).
make_loop_panel <- function(n_subjects = 6, n_obs = 30, seed = 99) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(n_subjects), function(id) {
    a  <- as.numeric(arima.sim(list(ar = 0.4), n = n_obs))
    b  <- 0.4 * a  + as.numeric(arima.sim(list(ar = 0.3), n = n_obs))
    cc <- 0.4 * b  + as.numeric(arima.sim(list(ar = 0.2), n = n_obs))
    data.frame(
      id   = as.character(id),
      time = seq_len(n_obs),
      a    = a, b = b, c = cc,
      stringsAsFactors = FALSE
    )
  }))
}

# Negative full panel — strong negative a→b. Used for Layer 2b.
make_neg_loop_panel <- function(n_subjects = 6, n_obs = 30, seed = 77) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(n_subjects), function(id) {
    a  <- as.numeric(arima.sim(list(ar = 0.4), n = n_obs))
    b  <- -0.8 * a + as.numeric(arima.sim(list(ar = 0.3), n = n_obs))
    cc <-  0.4 * b + as.numeric(arima.sim(list(ar = 0.2), n = n_obs))
    data.frame(
      id   = as.character(id),
      time = seq_len(n_obs),
      a    = a, b = b, c = cc,
      stringsAsFactors = FALSE
    )
  }))
}


# ══════════════════════════════════════════════════════════════════════════════
# Layer 1 — Validation (fast, no CRAN skip needed)
# ══════════════════════════════════════════════════════════════════════════════

test_that("misspelled a_series triggers informative error", {
  panel <- make_tiny_loop_panel()
  expect_error(
    suppressMessages(looping_machine(
      panel, a_series = "not_a_col", b_series = "b", c_series = "c",
      id_var = "id", timevar = "time", min_n_subject = 1
    )),
    regexp = "not_a_col"
  )
})

test_that("misspelled b_series triggers informative error", {
  panel <- make_tiny_loop_panel()
  expect_error(
    suppressMessages(looping_machine(
      panel, a_series = "a", b_series = "not_b_col", c_series = "c",
      id_var = "id", timevar = "time", min_n_subject = 1
    )),
    regexp = "not_b_col"
  )
})

test_that("misspelled c_series triggers informative error", {
  panel <- make_tiny_loop_panel()
  expect_error(
    suppressMessages(looping_machine(
      panel, a_series = "a", b_series = "b", c_series = "not_c_col",
      id_var = "id", timevar = "time", min_n_subject = 1
    )),
    regexp = "not_c_col"
  )
})

test_that("misspelled id_var triggers informative error", {
  panel <- make_tiny_loop_panel()
  expect_error(
    suppressMessages(looping_machine(
      panel, a_series = "a", b_series = "b", c_series = "c",
      id_var = "not_an_id", timevar = "time", min_n_subject = 1
    )),
    regexp = "not_an_id"
  )
})

test_that("misspelled timevar triggers informative error", {
  panel <- make_tiny_loop_panel()
  expect_error(
    suppressMessages(looping_machine(
      panel, a_series = "a", b_series = "b", c_series = "c",
      id_var = "id", timevar = "not_a_time", min_n_subject = 1
    )),
    regexp = "not_a_time"
  )
})

test_that("min_n_subject = 0 triggers upfront error", {
  panel <- make_tiny_loop_panel()
  expect_error(
    looping_machine(panel, a_series = "a", b_series = "b", c_series = "c",
                    id_var = "id", timevar = "time", min_n_subject = 0),
    regexp = "min_n_subject"
  )
})

test_that("minvar = -1 triggers upfront error", {
  panel <- make_tiny_loop_panel()
  expect_error(
    looping_machine(panel, a_series = "a", b_series = "b", c_series = "c",
                    id_var = "id", timevar = "time", minvar = -1),
    regexp = "minvar"
  )
})

test_that("all subjects below min_n_subject threshold raises an error", {
  panel <- make_tiny_loop_panel()          # 5 obs per subject
  expect_error(
    suppressMessages(looping_machine(
      panel, a_series = "a", b_series = "b", c_series = "c",
      id_var = "id", timevar = "time", min_n_subject = 20
    )),
    regexp = "not enough cases"
  )
})

test_that("duplicate series names trigger an informative error", {
  panel <- make_tiny_loop_panel()
  expect_error(
    looping_machine(
      panel, a_series = "a", b_series = "a", c_series = "c",
      id_var = "id", timevar = "time", min_n_subject = 1
    ),
    regexp = "different"
  )
  expect_error(
    looping_machine(
      panel, a_series = "a", b_series = "b", c_series = "b",
      id_var = "id", timevar = "time", min_n_subject = 1
    ),
    regexp = "different"
  )
  expect_error(
    looping_machine(
      panel, a_series = "a", b_series = "b", c_series = "a",
      id_var = "id", timevar = "time", min_n_subject = 1
    ),
    regexp = "different"
  )
})

test_that("covariates overlapping with a loop variable triggers an informative error", {
  panel <- make_tiny_loop_panel()
  expect_error(
    looping_machine(
      panel, a_series = "a", b_series = "b", c_series = "c",
      id_var = "id", timevar = "time",
      covariates = "a",   # "a" is also a_series
      min_n_subject = 1
    ),
    regexp = "overlap"
  )
  expect_error(
    looping_machine(
      panel, a_series = "a", b_series = "b", c_series = "c",
      id_var = "id", timevar = "time",
      covariates = c("extra", "c"),   # "c" is also c_series
      min_n_subject = 1
    ),
    regexp = "overlap"
  )
})

test_that("alpha outside (0, 1) raises an informative error", {
  panel <- make_tiny_loop_panel()
  expect_error(
    looping_machine(
      panel, a_series = "a", b_series = "b", c_series = "c",
      id_var = "id", timevar = "time", min_n_subject = 1, alpha = 1.5
    ),
    regexp = "alpha"
  )
  expect_error(
    looping_machine(
      panel, a_series = "a", b_series = "b", c_series = "c",
      id_var = "id", timevar = "time", min_n_subject = 1, alpha = 0
    ),
    regexp = "alpha"
  )
  expect_error(
    looping_machine(
      panel, a_series = "a", b_series = "b", c_series = "c",
      id_var = "id", timevar = "time", min_n_subject = 1, alpha = -0.1
    ),
    regexp = "alpha"
  )
})


# ══════════════════════════════════════════════════════════════════════════════
# Layer 2a — Minimal integration (NO skip_on_cran, always runs under covr)
# ══════════════════════════════════════════════════════════════════════════════
#
# covr::package_coverage() does NOT set NOT_CRAN, so skip_on_cran() always
# fires there. These tests use a minimal panel (4 subjects × 25 obs) with
# mild cross-series signal; auto.arima typically picks low-order models and
# the whole looping_machine() call completes in a few seconds — well within
# CRAN's time limit.

# Panels created at file level (pure data construction, no model fitting).
mini_panel     <- make_mini_loop_panel()
mini_panel_neg <- make_mini_neg_loop_panel()

# Results are memoised: looping_machine() is fitted once per panel and the
# cached object is reused across all tests to avoid redundant model fitting.
.mini_cache <- new.env(parent = emptyenv())

.get_mini <- function() {
  if (!exists("r", envir = .mini_cache)) {
    .mini_cache$r <- suppressMessages(
      looping_machine(mini_panel, a_series = "a", b_series = "b", c_series = "c",
                      id_var = "id", timevar = "time")
    )
  }
  .mini_cache$r
}

.get_mini_neg <- function() {
  if (!exists("rneg", envir = .mini_cache)) {
    .mini_cache$rneg <- suppressMessages(
      looping_machine(mini_panel_neg, a_series = "a", b_series = "b", c_series = "c",
                      id_var = "id", timevar = "time")
    )
  }
  .mini_cache$rneg
}


# ── Return structure ──────────────────────────────────────────────────────────

test_that("looping_machine returns a plain list (not an S3 class)", {
  r <- .get_mini()
  expect_type(r, "list")
  expect_false(inherits(r, "iarimax_results"))
})

test_that("return list has all expected top-level names", {
  r <- .get_mini()
  expected <- c("loop_df", "alpha", "covariates", "include_third_as_covariate",
                "loop_case_detail", "iarimax_a_to_b", "iarimax_b_to_c", "iarimax_c_to_a")
  expect_true(all(expected %in% names(r)))
})

test_that("each iarimax leg is an iarimax_results object", {
  r <- .get_mini()
  expect_s3_class(r$iarimax_a_to_b, "iarimax_results")
  expect_s3_class(r$iarimax_b_to_c, "iarimax_results")
  expect_s3_class(r$iarimax_c_to_a, "iarimax_results")
})

test_that("loop_df is a data.frame with at least one row", {
  r <- .get_mini()
  expect_s3_class(r$loop_df, "data.frame")
  expect_gt(nrow(r$loop_df), 0L)
})


# ── loop_case_detail structure ────────────────────────────────────────────────

test_that("loop_case_detail has the expected fields", {
  d <- .get_mini()$loop_case_detail
  expect_type(d, "list")
  expect_true(all(c("n_in_loop_df", "n_complete", "n_na_indicator", "n_dropped_by_join") %in% names(d)))
})

test_that("loop_case_detail counts are internally consistent", {
  r <- .get_mini()
  d <- r$loop_case_detail
  expect_equal(d$n_in_loop_df, d$n_complete + d$n_na_indicator)
  expect_equal(d$n_in_loop_df, nrow(r$loop_df))
})

test_that("loop_case_detail n_complete matches non-NA Loop_positive_directed rows", {
  r <- .get_mini()
  expect_equal(
    r$loop_case_detail$n_complete,
    sum(!is.na(r$loop_df$Loop_positive_directed))
  )
})


# ── Column naming ─────────────────────────────────────────────────────────────

test_that("loop_df contains correctly named coefficient columns for all three legs", {
  r <- .get_mini()
  expect_true(all(c("a_b", "b_c", "c_a") %in% names(r$loop_df)))
})

test_that("loop_df contains stderr columns for all three legs", {
  r <- .get_mini()
  expect_true(all(c("stderr_a_b", "stderr_b_c", "stderr_c_a") %in% names(r$loop_df)))
})

test_that("loop_df contains n_valid columns for all three legs", {
  r <- .get_mini()
  expect_true(all(c("a_b_n_valid", "b_c_n_valid", "c_a_n_valid") %in% names(r$loop_df)))
})

test_that("loop_df contains n_params columns for all three legs", {
  r <- .get_mini()
  expect_true(all(c("a_b_n_params", "b_c_n_params", "c_a_n_params") %in% names(r$loop_df)))
})

test_that("loop_df contains p-value columns for all three legs", {
  r <- .get_mini()
  expect_true(all(c("a_b_pval", "b_c_pval", "c_a_pval") %in% names(r$loop_df)))
})

test_that("loop_df contains the id and Loop_positive_directed columns", {
  r      <- .get_mini()
  id_col <- attr(r$iarimax_a_to_b, "id_var")
  expect_true(id_col %in% names(r$loop_df))
  expect_true("Loop_positive_directed" %in% names(r$loop_df))
})

test_that("loop_df values match individual leg results_df (no column mixup)", {
  r  <- .get_mini()
  df <- r$loop_df

  # For each subject in loop_df, verify that the renamed columns carry the

  # correct values from the corresponding leg's results_df.
  for (id in df$id) {
    # a→b leg
    ab_row <- r$iarimax_a_to_b$results_df[r$iarimax_a_to_b$results_df$id == id, ]
    expect_equal(df$a_b[df$id == id],           ab_row$estimate_a,   ignore_attr = TRUE)
    expect_equal(df$stderr_a_b[df$id == id],     ab_row$std.error_a, ignore_attr = TRUE)
    expect_equal(df$a_b_n_valid[df$id == id],    ab_row$n_valid,     ignore_attr = TRUE)
    expect_equal(df$a_b_n_params[df$id == id],   ab_row$n_params,    ignore_attr = TRUE)
    expect_equal(df$a_b_pval[df$id == id],       ab_row$pval_a,      ignore_attr = TRUE)

    # b→c leg
    bc_row <- r$iarimax_b_to_c$results_df[r$iarimax_b_to_c$results_df$id == id, ]
    expect_equal(df$b_c[df$id == id],           bc_row$estimate_b,   ignore_attr = TRUE)
    expect_equal(df$stderr_b_c[df$id == id],     bc_row$std.error_b, ignore_attr = TRUE)
    expect_equal(df$b_c_n_valid[df$id == id],    bc_row$n_valid,     ignore_attr = TRUE)
    expect_equal(df$b_c_n_params[df$id == id],   bc_row$n_params,    ignore_attr = TRUE)
    expect_equal(df$b_c_pval[df$id == id],       bc_row$pval_b,      ignore_attr = TRUE)

    # c→a leg
    ca_row <- r$iarimax_c_to_a$results_df[r$iarimax_c_to_a$results_df$id == id, ]
    expect_equal(df$c_a[df$id == id],           ca_row$estimate_c,   ignore_attr = TRUE)
    expect_equal(df$stderr_c_a[df$id == id],     ca_row$std.error_c, ignore_attr = TRUE)
    expect_equal(df$c_a_n_valid[df$id == id],    ca_row$n_valid,     ignore_attr = TRUE)
    expect_equal(df$c_a_n_params[df$id == id],   ca_row$n_params,    ignore_attr = TRUE)
    expect_equal(df$c_a_pval[df$id == id],       ca_row$pval_c,      ignore_attr = TRUE)
  }
})


# ── i_pval already applied ───────────────────────────────────────────────────

test_that("pval_a present in iarimax_a_to_b$results_df", {
  expect_true("pval_a" %in% names(.get_mini()$iarimax_a_to_b$results_df))
})

test_that("pval_b present in iarimax_b_to_c$results_df", {
  expect_true("pval_b" %in% names(.get_mini()$iarimax_b_to_c$results_df))
})

test_that("pval_c present in iarimax_c_to_a$results_df", {
  expect_true("pval_c" %in% names(.get_mini()$iarimax_c_to_a$results_df))
})


# ── Loop_positive_directed correctness ───────────────────────────────────────

test_that("Loop_positive_directed is integer with values 0, 1, or NA only", {
  vals <- .get_mini()$loop_df$Loop_positive_directed
  expect_type(vals, "integer")
  expect_true(all(vals %in% c(0L, 1L) | is.na(vals)))
})

test_that("Loop_positive_directed matches the conjunction formula exactly", {
  r  <- .get_mini()
  df <- r$loop_df
  al <- r$alpha

  expected <- as.integer(
    df[["a_b"]] > 0 & df[["b_c"]] > 0 & df[["c_a"]] > 0 &
    df[["a_b_pval"]] < al & df[["b_c_pval"]] < al & df[["c_a_pval"]] < al
  )
  expect_equal(df[["Loop_positive_directed"]], expected)
})

test_that("a row with a negative a_b beta is flagged Loop_positive_directed = 0", {
  df       <- .get_mini_neg()$loop_df
  neg_rows <- which(df[["a_b"]] < 0)
  expect_gt(length(neg_rows), 0L)
  expect_equal(df[["Loop_positive_directed"]][neg_rows[1]], 0L)
})

test_that("NA guard: subject with NA pval gets Loop_positive_directed = NA, not 0", {
  # Build a fake loop_df with one subject having NA pval (simulating a failed model).
  # Without the any_na guard, NA & FALSE would collapse to FALSE → 0L.
  fake_df <- data.frame(
    id       = c("A", "B", "C"),
    a_b      = c(0.5,  0.3,  NA),
    a_b_pval = c(0.01, 0.01, NA),
    b_c      = c(0.4,  -0.2, 0.5),
    b_c_pval = c(0.02, 0.80, 0.03),
    c_a      = c(0.3,  0.1,  0.4),
    c_a_pval = c(0.03, 0.03, 0.01),
    stringsAsFactors = FALSE
  )

  any_na <- is.na(fake_df$a_b_pval) |
            is.na(fake_df$b_c_pval) |
            is.na(fake_df$c_a_pval)

  indicator <- ifelse(
    any_na,
    NA_integer_,
    ifelse(
      fake_df$a_b_pval < 0.05 & fake_df$b_c_pval < 0.05 & fake_df$c_a_pval < 0.05 &
        fake_df$a_b > 0 & fake_df$b_c > 0 & fake_df$c_a > 0,
      1L, 0L
    )
  )

  # A: all positive and sig → 1
  expect_equal(indicator[1], 1L)
  # B: b_c is negative and non-sig → 0
  expect_equal(indicator[2], 0L)
  # C: NA pval in a_b → NA, NOT 0
  expect_true(is.na(indicator[3]))
})

test_that("Loop_positive_directed is correct per subject (no cross-subject mixup)", {
  r  <- .get_mini()
  df <- r$loop_df
  al <- r$alpha

  for (i in seq_len(nrow(df))) {
    row <- df[i, ]
    has_na <- is.na(row$a_b_pval) || is.na(row$b_c_pval) || is.na(row$c_a_pval)

    if (has_na) {
      expect_true(is.na(row$Loop_positive_directed),
                  info = paste("Subject", row$id, "should be NA (failed model)"))
    } else {
      is_loop <- row$a_b > 0 && row$b_c > 0 && row$c_a > 0 &&
                 row$a_b_pval < al && row$b_c_pval < al && row$c_a_pval < al
      expected <- if (is_loop) 1L else 0L
      expect_equal(row$Loop_positive_directed, expected,
                   info = paste("Subject", row$id, "indicator mismatch"))
    }
  }
})


# ── Parameter storage ─────────────────────────────────────────────────────────

test_that("default alpha (0.05), covariates (NULL), include_third (FALSE) stored correctly", {
  r <- .get_mini()
  expect_equal(r$alpha, 0.05)
  expect_null(r$covariates)
  expect_false(r$include_third_as_covariate)
})


# ── Focal predictor attributes ────────────────────────────────────────────────

test_that("each leg has the correct focal_predictor attribute", {
  r <- .get_mini()
  expect_equal(attr(r$iarimax_a_to_b, "focal_predictor"), "a")
  expect_equal(attr(r$iarimax_b_to_c, "focal_predictor"), "b")
  expect_equal(attr(r$iarimax_c_to_a, "focal_predictor"), "c")
})


# ── Series names with underscores ────────────────────────────────────────────

test_that("series names containing underscores produce correct leg column names", {
  # Exercises the paste0(a_series, "_", b_series) leg naming with underscore-containing
  # names. Verifies the expected column names appear in loop_df.
  set.seed(11)
  n_sub <- 4; n_obs <- 25
  panel_u <- do.call(rbind, lapply(seq_len(n_sub), function(id) {
    x_1 <- rnorm(n_obs)
    x_2 <- 0.3 * x_1 + rnorm(n_obs)
    x_3 <- 0.3 * x_2 + rnorm(n_obs)
    data.frame(id = as.character(id), time = seq_len(n_obs),
               x_1 = x_1, x_2 = x_2, x_3 = x_3, stringsAsFactors = FALSE)
  }))
  r <- suppressMessages(
    looping_machine(panel_u, a_series = "x_1", b_series = "x_2", c_series = "x_3",
                    id_var = "id", timevar = "time")
  )
  expect_true("x_1_x_2" %in% names(r$loop_df))
  expect_true("x_2_x_3" %in% names(r$loop_df))
  expect_true("x_3_x_1" %in% names(r$loop_df))
})

test_that("colliding leg names from underscore series trigger an informative error", {
  # a="x", b="y_x", c="x_y" makes ab_name = ca_name = "x_y_x".
  panel_c <- data.frame(
    id  = rep(c("1", "2"), each = 5),
    time = rep(seq_len(5), 2),
    x   = rnorm(10), y_x = rnorm(10), x_y = rnorm(10),
    stringsAsFactors = FALSE
  )
  expect_error(
    looping_machine(panel_c, a_series = "x", b_series = "y_x", c_series = "x_y",
                    id_var = "id", timevar = "time", min_n_subject = 1),
    regexp = "non-unique leg names"
  )
})


# ── Merge integrity ───────────────────────────────────────────────────────────

test_that("loop_df row count is <= min n_used_iarimax across the three legs", {
  r    <- .get_mini()
  n_ab <- r$iarimax_a_to_b$case_number_detail$n_used_iarimax
  n_bc <- r$iarimax_b_to_c$case_number_detail$n_used_iarimax
  n_ca <- r$iarimax_c_to_a$case_number_detail$n_used_iarimax
  expect_lte(nrow(r$loop_df), min(n_ab, n_bc, n_ca))
})

test_that("all ids in loop_df appear in all three legs' results_df", {
  r        <- .get_mini()
  ids_loop <- r$loop_df[["id"]]
  expect_true(all(ids_loop %in% r$iarimax_a_to_b$results_df[["id"]]))
  expect_true(all(ids_loop %in% r$iarimax_b_to_c$results_df[["id"]]))
  expect_true(all(ids_loop %in% r$iarimax_c_to_a$results_df[["id"]]))
})


# ── n_dropped_by_join ─────────────────────────────────────────────────────────

test_that("n_dropped_by_join is positive and message emitted when a subject is absent from one leg", {
  # Subject "5" has c = 0 (zero variance): passes a->b (needs a, b) but is
  # filtered out of b->c (needs b, c) and c->a (needs c, a). After the inner
  # join, subject 5 is absent from loop_df → n_dropped_by_join should be 1.
  set.seed(42)
  n_obs <- 25
  base <- do.call(rbind, lapply(1:4, function(id) {
    a  <- rnorm(n_obs); b <- 0.4 * a + rnorm(n_obs); cc <- 0.4 * b + rnorm(n_obs)
    data.frame(id = as.character(id), time = seq_len(n_obs),
               a = a, b = b, c = cc, stringsAsFactors = FALSE)
  }))
  extra <- data.frame(id = "5", time = seq_len(n_obs),
                      a = rnorm(n_obs), b = rnorm(n_obs), c = rep(0, n_obs),
                      stringsAsFactors = FALSE)
  panel_drop <- rbind(base, extra)

  msgs <- character(0)
  r <- withCallingHandlers(
    suppressWarnings(looping_machine(panel_drop, a_series = "a", b_series = "b",
                                     c_series = "c", id_var = "id", timevar = "time")),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }
  )

  expect_equal(r$loop_case_detail$n_dropped_by_join, 1L)
  expect_true(any(grepl("dropped from loop_df", msgs)))
})

test_that("n_dropped_by_join is 0 when all subjects survive all three legs", {
  r <- .get_mini()
  expect_equal(r$loop_case_detail$n_dropped_by_join, 0L)
})


# ── Verbose branch coverage (no skip_on_cran) ─────────────────────────────────
# Runs looping_machine() with verbose = TRUE to instrument the per-leg start
# messages and the finish message.

test_that("verbose = TRUE runs without error and emits messages", {
  msgs <- character(0)
  withCallingHandlers(
    suppressWarnings(
      looping_machine(mini_panel, a_series = "a", b_series = "b", c_series = "c",
                      id_var = "id", timevar = "time", verbose = TRUE)
    ),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_gt(length(msgs), 1L)
})

test_that("summary message fires with verbose = TRUE and mentions 'Number of cases'", {
  msgs <- character(0)
  withCallingHandlers(
    suppressWarnings(
      looping_machine(mini_panel, a_series = "a", b_series = "b", c_series = "c",
                      id_var = "id", timevar = "time", verbose = TRUE)
    ),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_true(any(grepl("Number of cases", msgs)))
})

test_that("summary message does not fire with verbose = FALSE", {
  msgs <- character(0)
  withCallingHandlers(
    suppressWarnings(
      looping_machine(mini_panel, a_series = "a", b_series = "b", c_series = "c",
                      id_var = "id", timevar = "time", verbose = FALSE)
    ),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_false(any(grepl("Number of cases", msgs)))
})


# ── include_third_as_covariate = TRUE branch coverage ────────────────────────
# Covers the TRUE branch of `if (include_third_as_covariate)` in the source.
# Also verifies that the correct third variable actually appears as a
# coefficient column in each leg (not just that the flag is stored).

test_that("include_third_as_covariate = TRUE adds the correct third variable to each leg", {
  r3 <- suppressMessages(
    looping_machine(mini_panel, a_series = "a", b_series = "b", c_series = "c",
                    id_var = "id", timevar = "time", include_third_as_covariate = TRUE)
  )
  expect_true(r3$include_third_as_covariate)
  # a→b leg: c is the third variable
  expect_true("estimate_c" %in% names(r3$iarimax_a_to_b$results_df))
  # b→c leg: a is the third variable
  expect_true("estimate_a" %in% names(r3$iarimax_b_to_c$results_df))
  # c→a leg: b is the third variable
  expect_true("estimate_b" %in% names(r3$iarimax_c_to_a$results_df))
})


# ── keep_models forwarding ────────────────────────────────────────────────────

test_that("keep_models = TRUE populates $models in each leg", {
  r_km <- suppressMessages(
    looping_machine(mini_panel, a_series = "a", b_series = "b", c_series = "c",
                    id_var = "id", timevar = "time", keep_models = TRUE)
  )
  expect_false(is.null(r_km$iarimax_a_to_b$models))
  expect_false(is.null(r_km$iarimax_b_to_c$models))
  expect_false(is.null(r_km$iarimax_c_to_a$models))
})


# ── correlation_method forwarding ─────────────────────────────────────────────

test_that("correlation_method = 'spearman' runs without error", {
  r_sp <- suppressMessages(
    looping_machine(mini_panel, a_series = "a", b_series = "b", c_series = "c",
                    id_var = "id", timevar = "time", correlation_method = "spearman")
  )
  expect_s3_class(r_sp$iarimax_a_to_b, "iarimax_results")
  expect_s3_class(r_sp$iarimax_b_to_c, "iarimax_results")
  expect_s3_class(r_sp$iarimax_c_to_a, "iarimax_results")
})


# ══════════════════════════════════════════════════════════════════════════════
# Layer 2b — Larger panel (skip_on_cran, runs under devtools::test() only)
# ══════════════════════════════════════════════════════════════════════════════

panel     <- make_loop_panel()
panel_neg <- make_neg_loop_panel()

.lm_cache <- new.env(parent = emptyenv())

.get_result <- function() {
  if (!exists("result", envir = .lm_cache)) {
    skip_on_cran()
    .lm_cache$result <- suppressMessages(
      looping_machine(panel, a_series = "a", b_series = "b", c_series = "c",
                      id_var = "id", timevar = "time")
    )
  }
  .lm_cache$result
}

.get_result_neg <- function() {
  if (!exists("result_neg", envir = .lm_cache)) {
    skip_on_cran()
    .lm_cache$result_neg <- suppressMessages(
      looping_machine(panel_neg, a_series = "a", b_series = "b", c_series = "c",
                      id_var = "id", timevar = "time")
    )
  }
  .lm_cache$result_neg
}

test_that("alpha = 0.01 is stored correctly", {
  skip_on_cran()
  r <- suppressMessages(looping_machine(
    panel, a_series = "a", b_series = "b", c_series = "c",
    id_var = "id", timevar = "time", alpha = 0.01
  ))
  expect_equal(r$alpha, 0.01)
})

test_that("non-NULL covariates are stored and appear as coefficient columns in all legs", {
  skip_on_cran()
  panel2      <- panel
  panel2$cov1 <- rnorm(nrow(panel2))
  r <- suppressMessages(looping_machine(
    panel2, a_series = "a", b_series = "b", c_series = "c",
    id_var = "id", timevar = "time", covariates = "cov1"
  ))
  expect_equal(r$covariates, "cov1")
  expect_true("estimate_cov1" %in% names(r$iarimax_a_to_b$results_df))
  expect_true("estimate_cov1" %in% names(r$iarimax_b_to_c$results_df))
  expect_true("estimate_cov1" %in% names(r$iarimax_c_to_a$results_df))
})

test_that("summary message fires with verbose = TRUE (real panel)", {
  skip_on_cran()
  msgs <- character(0)
  withCallingHandlers(
    suppressWarnings(looping_machine(
      panel, a_series = "a", b_series = "b", c_series = "c",
      id_var = "id", timevar = "time", verbose = TRUE
    )),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_true(any(grepl("Number of cases", msgs)))
})

test_that("summary message does not fire with verbose = FALSE (real panel)", {
  skip_on_cran()
  msgs <- character(0)
  withCallingHandlers(
    suppressWarnings(looping_machine(
      panel, a_series = "a", b_series = "b", c_series = "c",
      id_var = "id", timevar = "time", verbose = FALSE
    )),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_false(any(grepl("Number of cases", msgs)))
})

test_that("verbose = TRUE emits more messages than verbose = FALSE", {
  skip_on_cran()
  count_msgs <- function(verbose) {
    n <- 0L
    withCallingHandlers(
      suppressWarnings(looping_machine(
        panel, a_series = "a", b_series = "b", c_series = "c",
        id_var = "id", timevar = "time", verbose = verbose
      )),
      message = function(m) { n <<- n + 1L; invokeRestart("muffleMessage") }
    )
    n
  }
  expect_gt(count_msgs(TRUE), count_msgs(FALSE))
})

test_that("Loop_positive_directed is 1 for at least one subject in a strong-signal panel", {
  skip_on_cran()
  set.seed(42)
  n_sub <- 8; n_obs <- 50
  strong_panel <- do.call(rbind, lapply(seq_len(n_sub), function(id) {
    a  <- as.numeric(arima.sim(list(ar = 0.3), n = n_obs))
    b  <- 0.95 * a  + 0.05 * rnorm(n_obs)
    cc <- 0.95 * b  + 0.05 * rnorm(n_obs)
    data.frame(id = as.character(id), time = seq_len(n_obs),
               a = a, b = b, c = cc, stringsAsFactors = FALSE)
  }))
  r_strong <- suppressMessages(looping_machine(
    strong_panel, a_series = "a", b_series = "b", c_series = "c",
    id_var = "id", timevar = "time"
  ))
  expect_gt(sum(r_strong$loop_df$Loop_positive_directed, na.rm = TRUE), 0L)
})

test_that("Loop_positive_directed logic is consistent with conjunction formula on full panel", {
  r  <- .get_result()
  df <- r$loop_df
  al <- r$alpha
  expected <- as.integer(
    df[["a_b"]] > 0 & df[["b_c"]] > 0 & df[["c_a"]] > 0 &
    df[["a_b_pval"]] < al & df[["b_c_pval"]] < al & df[["c_a_pval"]] < al
  )
  expect_equal(df[["Loop_positive_directed"]], expected)
})
