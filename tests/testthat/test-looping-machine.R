# Tests for looping_machine()
# Layer 1: validation and structural checks using a tiny dataset (errors thrown
#           before auto.arima runs) — no skip_on_cran() needed.
# Layer 2: full integration tests on a real panel — skip_on_cran().

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

# Full panel — 6 subjects × 30 obs with mild AR(1) dynamics and cross-series
# signal so auto.arima reliably selects low-order models and the loop indicator
# can in principle fire.
make_loop_panel <- function(n_subjects = 6, n_obs = 30, seed = 99) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(n_subjects), function(id) {
    a <- as.numeric(arima.sim(list(ar = 0.4), n = n_obs))
    b <- 0.4 * a + as.numeric(arima.sim(list(ar = 0.3), n = n_obs))
    cc <- 0.4 * b + as.numeric(arima.sim(list(ar = 0.2), n = n_obs))
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

# Strong-signal panel — 20 subjects × 80 obs with b = 0.9*a and c = 0.9*b,
# so c ≈ 0.81*a. All three legs (a→b, b→c, c→a) should reliably show
# significant positive coefficients, ensuring Loop_positive_directed = 1
# fires for at least some subjects.
make_strong_loop_panel <- function(n_subjects = 20, n_obs = 80, seed = 42) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(n_subjects), function(id) {
    a  <- as.numeric(arima.sim(list(ar = 0.3), n = n_obs))
    b  <- 0.9 * a + 0.2 * rnorm(n_obs)
    cc <- 0.9 * b + 0.2 * rnorm(n_obs)
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
# Layer 2 — Integration tests (real model fitting, skip on CRAN)
# ══════════════════════════════════════════════════════════════════════════════

skip_on_cran()

panel <- make_loop_panel()
result <- suppressMessages(looping_machine(
  panel, a_series = "a", b_series = "b", c_series = "c",
  id_var = "id", timevar = "time"
))

# Panel with a strong negative a→b relationship so at least one subject's
# estimated a_b coefficient is reliably negative.
make_neg_loop_panel <- function(n_subjects = 6, n_obs = 30, seed = 77) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(n_subjects), function(id) {
    a  <- as.numeric(arima.sim(list(ar = 0.4), n = n_obs))
    b  <- -0.8 * a + as.numeric(arima.sim(list(ar = 0.3), n = n_obs))
    cc <- 0.4 * b  + as.numeric(arima.sim(list(ar = 0.2), n = n_obs))
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
panel_neg <- make_neg_loop_panel()
result_neg <- suppressMessages(looping_machine(
  panel_neg, a_series = "a", b_series = "b", c_series = "c",
  id_var = "id", timevar = "time"
))

# Strong-signal panel: used to verify Loop_positive_directed = 1 can actually fire.
panel_strong <- make_strong_loop_panel()
result_strong <- suppressMessages(looping_machine(
  panel_strong, a_series = "a", b_series = "b", c_series = "c",
  id_var = "id", timevar = "time"
))


# ── Return structure ──────────────────────────────────────────────────────────

test_that("looping_machine returns a plain list (not an S3 class)", {
  expect_type(result, "list")
  expect_false(inherits(result, "iarimax_results"))
})

test_that("return list has all expected top-level names", {
  expected_names <- c(
    "loop_df", "alpha", "covariates", "include_third_as_covariate",
    "loop_case_detail",
    "iarimax_a_to_b", "iarimax_b_to_c", "iarimax_c_to_a"
  )
  expect_true(all(expected_names %in% names(result)))
})

test_that("each iarimax leg is an iarimax_results object", {
  expect_s3_class(result$iarimax_a_to_b, "iarimax_results")
  expect_s3_class(result$iarimax_b_to_c, "iarimax_results")
  expect_s3_class(result$iarimax_c_to_a, "iarimax_results")
})

test_that("loop_df is a data.frame", {
  expect_s3_class(result$loop_df, "data.frame")
})

test_that("loop_df has at least one row", {
  expect_gt(nrow(result$loop_df), 0L)
})


# ── loop_case_detail structure ────────────────────────────────────────────────

test_that("loop_case_detail is a list with the expected fields", {
  d <- result$loop_case_detail
  expect_type(d, "list")
  expect_true(all(c("n_in_loop_df", "n_complete", "n_na_indicator") %in% names(d)))
})

test_that("loop_case_detail counts are internally consistent", {
  d <- result$loop_case_detail
  expect_equal(d$n_in_loop_df, d$n_complete + d$n_na_indicator)
  expect_equal(d$n_in_loop_df, nrow(result$loop_df))
})

test_that("loop_case_detail n_complete matches non-NA Loop_positive_directed rows", {
  expect_equal(
    result$loop_case_detail$n_complete,
    sum(!is.na(result$loop_df$Loop_positive_directed))
  )
})


# ── Column naming ─────────────────────────────────────────────────────────────

test_that("loop_df contains correctly named coefficient columns for all three legs", {
  expected_cols <- c("a_b", "b_c", "c_a")
  expect_true(all(expected_cols %in% names(result$loop_df)))
})

test_that("loop_df contains stderr columns for all three legs", {
  expected_cols <- c("stderr_a_b", "stderr_b_c", "stderr_c_a")
  expect_true(all(expected_cols %in% names(result$loop_df)))
})

test_that("loop_df contains n_valid columns for all three legs", {
  expected_cols <- c("a_b_n_valid", "b_c_n_valid", "c_a_n_valid")
  expect_true(all(expected_cols %in% names(result$loop_df)))
})

test_that("loop_df contains p-value columns for all three legs", {
  expected_cols <- c("a_b_pval", "b_c_pval", "c_a_pval")
  expect_true(all(expected_cols %in% names(result$loop_df)))
})

test_that("loop_df contains the id column", {
  expect_true("id" %in% names(result$loop_df))
})

test_that("loop_df contains Loop_positive_directed column", {
  expect_true("Loop_positive_directed" %in% names(result$loop_df))
})


# ── i_pval already applied ───────────────────────────────────────────────────

test_that("pval column for a is present in iarimax_a_to_b$results_df", {
  expect_true("pval_a" %in% names(result$iarimax_a_to_b$results_df))
})

test_that("pval column for b is present in iarimax_b_to_c$results_df", {
  expect_true("pval_b" %in% names(result$iarimax_b_to_c$results_df))
})

test_that("pval column for c is present in iarimax_c_to_a$results_df", {
  expect_true("pval_c" %in% names(result$iarimax_c_to_a$results_df))
})


# ── Loop_positive_directed correctness ───────────────────────────────────────

test_that("Loop_positive_directed is integer", {
  expect_type(result$loop_df$Loop_positive_directed, "integer")
})

test_that("Loop_positive_directed only takes values 0, 1, or NA", {
  vals <- result$loop_df$Loop_positive_directed
  expect_true(all(vals %in% c(0L, 1L) | is.na(vals)))
})

test_that("Loop_positive_directed is 1 only when all three betas are positive and all three pvals < alpha", {
  df    <- result$loop_df
  alpha <- result$alpha

  is_positive_loop <- df[["a_b"]] > 0 & df[["b_c"]] > 0 & df[["c_a"]] > 0 &
    df[["a_b_pval"]] < alpha & df[["b_c_pval"]] < alpha & df[["c_a_pval"]] < alpha

  expect_equal(
    df[["Loop_positive_directed"]],
    as.integer(is_positive_loop)
  )
})

test_that("a row with any negative focal beta is flagged as Loop_positive_directed = 0", {
  # result_neg uses b = -0.8*a + noise, so a_b betas are reliably negative.
  df       <- result_neg$loop_df
  neg_rows <- which(df[["a_b"]] < 0)
  expect_gt(length(neg_rows), 0L)   # at least one negative beta must exist
  expect_equal(df[["Loop_positive_directed"]][neg_rows[1]], 0L)
})

test_that("Loop_positive_directed is 1 for at least one subject in the strong-signal panel", {
  # panel_strong has b = 0.9*a + noise and c = 0.9*b + noise, so c ≈ 0.81*a.
  # All three legs should show significant positive effects for most subjects.
  expect_gt(
    sum(result_strong$loop_df$Loop_positive_directed, na.rm = TRUE),
    0L
  )
})


# ── Summary message always fires ──────────────────────────────────────────────

test_that("summary message is always emitted regardless of verbose", {
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
  expect_true(any(grepl("Number of cases", msgs)))
})


# ── Parameter storage ─────────────────────────────────────────────────────────

test_that("alpha parameter is stored correctly in output", {
  result_alpha <- suppressMessages(looping_machine(
    panel, a_series = "a", b_series = "b", c_series = "c",
    id_var = "id", timevar = "time", alpha = 0.01
  ))
  expect_equal(result_alpha$alpha, 0.01)
})

test_that("covariates = NULL is stored as NULL", {
  expect_null(result$covariates)
})

test_that("include_third_as_covariate = FALSE is stored correctly", {
  expect_false(result$include_third_as_covariate)
})

test_that("non-NULL covariates vector is stored correctly", {
  panel2 <- panel
  panel2$cov1 <- rnorm(nrow(panel2))
  result2 <- suppressMessages(looping_machine(
    panel2, a_series = "a", b_series = "b", c_series = "c",
    id_var = "id", timevar = "time", covariates = "cov1"
  ))
  expect_equal(result2$covariates, "cov1")
})

test_that("include_third_as_covariate = TRUE is stored correctly", {
  result3 <- suppressMessages(looping_machine(
    panel, a_series = "a", b_series = "b", c_series = "c",
    id_var = "id", timevar = "time", include_third_as_covariate = TRUE
  ))
  expect_true(result3$include_third_as_covariate)
})


# ── Include third as covariate smoke test ─────────────────────────────────────

test_that("include_third_as_covariate = TRUE runs without error and returns correct structure", {
  result3 <- suppressMessages(looping_machine(
    panel, a_series = "a", b_series = "b", c_series = "c",
    id_var = "id", timevar = "time", include_third_as_covariate = TRUE
  ))
  expect_s3_class(result3$iarimax_a_to_b, "iarimax_results")
  expect_true("Loop_positive_directed" %in% names(result3$loop_df))
})


# ── Merge integrity ───────────────────────────────────────────────────────────

test_that("loop_df row count is <= minimum n_used_iarimax across the three legs", {
  n_ab <- result$iarimax_a_to_b$case_number_detail$n_used_iarimax
  n_bc <- result$iarimax_b_to_c$case_number_detail$n_used_iarimax
  n_ca <- result$iarimax_c_to_a$case_number_detail$n_used_iarimax
  expect_lte(nrow(result$loop_df), min(n_ab, n_bc, n_ca))
})

test_that("all ids in loop_df appear in all three legs' results_df", {
  ids_loop <- result$loop_df[["id"]]
  ids_ab   <- result$iarimax_a_to_b$results_df[["id"]]
  ids_bc   <- result$iarimax_b_to_c$results_df[["id"]]
  ids_ca   <- result$iarimax_c_to_a$results_df[["id"]]
  expect_true(all(ids_loop %in% ids_ab))
  expect_true(all(ids_loop %in% ids_bc))
  expect_true(all(ids_loop %in% ids_ca))
})


# ── Focal predictor attributes ────────────────────────────────────────────────

test_that("iarimax_a_to_b has focal_predictor attribute set to a_series", {
  expect_equal(attr(result$iarimax_a_to_b, "focal_predictor"), "a")
})

test_that("iarimax_b_to_c has focal_predictor attribute set to b_series", {
  expect_equal(attr(result$iarimax_b_to_c, "focal_predictor"), "b")
})

test_that("iarimax_c_to_a has focal_predictor attribute set to c_series", {
  expect_equal(attr(result$iarimax_c_to_a, "focal_predictor"), "c")
})


# ── Verbose ───────────────────────────────────────────────────────────────────

test_that("verbose = TRUE emits more messages than verbose = FALSE", {
  count_msgs <- function(verbose) {
    msgs <- character(0)
    withCallingHandlers(
      suppressWarnings(looping_machine(
        panel, a_series = "a", b_series = "b", c_series = "c",
        id_var = "id", timevar = "time", verbose = verbose
      )),
      message = function(m) {
        msgs <<- c(msgs, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    )
    length(msgs)
  }
  expect_gt(count_msgs(TRUE), count_msgs(FALSE))
})

test_that("verbose = FALSE suppresses progress messages (only the summary message fires)", {
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
  expect_equal(length(msgs), 1L)
  expect_match(msgs[1], "Number of cases")
})
