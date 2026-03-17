# Fast tests: all errors are thrown before auto.arima is ever invoked.
# No skip_on_cran() needed — no model fitting happens here.

panel <- make_panel(n_subjects = 2)

# ── Missing column names ──────────────────────────────────────────────────────

test_that("misspelled y_series triggers informative error", {
  expect_error(
    iarimax(dataframe = panel, y_series = "not_a_col", x_series = "x",
            id_var = "id", timevar = "time"),
    regexp = "not_a_col"
  )
})

test_that("misspelled x_series triggers informative error", {
  expect_error(
    iarimax(dataframe = panel, y_series = "y", x_series = "not_a_col",
            id_var = "id", timevar = "time"),
    regexp = "not_a_col"
  )
})

test_that("misspelled id_var triggers informative error", {
  expect_error(
    iarimax(dataframe = panel, y_series = "y", x_series = "x",
            id_var = "not_an_id", timevar = "time"),
    regexp = "not_an_id"
  )
})

test_that("misspelled timevar triggers informative error", {
  expect_error(
    iarimax(dataframe = panel, y_series = "y", x_series = "x",
            id_var = "id", timevar = "not_a_time"),
    regexp = "not_a_time"
  )
})

test_that("error message includes 'Cannot find required variables'", {
  expect_error(
    iarimax(dataframe = panel, y_series = "oops", x_series = "x",
            id_var = "id", timevar = "time"),
    regexp = "Cannot find required variables"
  )
})

# ── NA in timevar ─────────────────────────────────────────────────────────────

test_that("single NA in timevar triggers error mentioning the count", {
  bad <- panel
  bad$time[1] <- NA
  expect_error(
    iarimax(dataframe = bad, y_series = "y", x_series = "x",
            id_var = "id", timevar = "time"),
    regexp = "1 row"
  )
})

test_that("single NA in timevar error names the offending variable", {
  bad <- panel
  bad$time[1] <- NA
  expect_error(
    iarimax(dataframe = bad, y_series = "y", x_series = "x",
            id_var = "id", timevar = "time"),
    regexp = "'time'"
  )
})

test_that("multiple NAs in timevar: error reports the exact count", {
  bad <- panel
  bad$time[c(1, 3, 7)] <- NA
  expect_error(
    iarimax(dataframe = bad, y_series = "y", x_series = "x",
            id_var = "id", timevar = "time"),
    regexp = "^3"
  )
})

test_that("NA in y does not trigger the timevar validation error", {
  bad <- panel
  bad$y[1] <- NA
  err <- tryCatch(
    iarimax(dataframe = bad, y_series = "y", x_series = "x",
            id_var = "id", timevar = "time"),
    error = function(e) e
  )
  if (inherits(err, "error")) {
    expect_false(
      grepl("missing values in the time variable", conditionMessage(err)),
      label = "NA in y must not trigger the timevar guard"
    )
  } else {
    succeed("iarimax ran successfully with NA in y")
  }
})

# ── focal_predictor validation ────────────────────────────────────────────────

test_that("multiple x_series without focal_predictor triggers error", {
  panel2 <- panel
  panel2$x2 <- rnorm(nrow(panel))
  expect_error(
    iarimax(dataframe = panel2, y_series = "y", x_series = c("x", "x2"),
            id_var = "id", timevar = "time"),
    regexp = "focal_predictor is required"
  )
})

test_that("focal_predictor not in x_series triggers error", {
  panel2 <- panel
  panel2$x2 <- rnorm(nrow(panel))
  expect_error(
    iarimax(dataframe = panel2, y_series = "y", x_series = c("x", "x2"),
            focal_predictor = "oops", id_var = "id", timevar = "time"),
    regexp = "focal_predictor must be one of"
  )
})

# ── correlation_method validation ────────────────────────────────────────────

test_that("invalid correlation_method triggers upfront error before the loop", {
  expect_error(
    iarimax(dataframe = panel, y_series = "y", x_series = "x",
            id_var = "id", timevar = "time", correlation_method = "blah"),
    regexp = "Correlation method not supported"
  )
})

test_that("invalid correlation_method error names the offending value", {
  expect_error(
    iarimax(dataframe = panel, y_series = "y", x_series = "x",
            id_var = "id", timevar = "time", correlation_method = "blah"),
    regexp = "blah"
  )
})

# ── all subjects filtered before loop ─────────────────────────────────────────

test_that("all subjects below min_n_subject threshold raises an error", {
  tiny <- data.frame(
    id   = rep(c("a", "b"), each = 5),
    time = rep(seq_len(5), 2),
    x    = rnorm(10),
    y    = rnorm(10),
    stringsAsFactors = FALSE
  )
  expect_error(
    iarimax(dataframe = tiny, y_series = "y", x_series = "x",
            id_var = "id", timevar = "time", min_n_subject = 20),
    regexp = "not enough cases"
  )
})
