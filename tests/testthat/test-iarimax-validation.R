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

# ── Non-numeric timevar ──────────────────────────────────────────────────────

test_that("character timevar triggers informative error", {
  bad <- panel
  bad$time <- as.character(bad$time)
  expect_error(
    iarimax(dataframe = bad, y_series = "y", x_series = "x",
            id_var = "id", timevar = "time"),
    regexp = "must be numeric"
  )
})

test_that("non-numeric timevar error reports the actual class", {
  bad <- panel
  bad$time <- as.Date("2020-01-01") + bad$time
  expect_error(
    iarimax(dataframe = bad, y_series = "y", x_series = "x",
            id_var = "id", timevar = "time"),
    regexp = "Date"
  )
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

# ── min_n_subject / minvar parameter validation ───────────────────────────────

test_that("min_n_subject = 0 triggers upfront error", {
  expect_error(
    iarimax(dataframe = panel, y_series = "y", x_series = "x",
            id_var = "id", timevar = "time", min_n_subject = 0),
    regexp = "min_n_subject"
  )
})

test_that("min_n_subject = Inf triggers upfront error", {
  expect_error(
    iarimax(dataframe = panel, y_series = "y", x_series = "x",
            id_var = "id", timevar = "time", min_n_subject = Inf),
    regexp = "min_n_subject"
  )
})

test_that("minvar = -1 triggers upfront error", {
  expect_error(
    iarimax(dataframe = panel, y_series = "y", x_series = "x",
            id_var = "id", timevar = "time", minvar = -1),
    regexp = "minvar"
  )
})

test_that("minvar = Inf triggers upfront error", {
  expect_error(
    iarimax(dataframe = panel, y_series = "y", x_series = "x",
            id_var = "id", timevar = "time", minvar = Inf),
    regexp = "minvar"
  )
})

# ── fixed_d validation ───────────────────────────────────────────────────────

test_that("fixed_d = -1 triggers upfront error", {
  expect_error(
    iarimax(dataframe = panel, y_series = "y", x_series = "x",
            id_var = "id", timevar = "time", fixed_d = -1),
    regexp = "fixed_d"
  )
})

test_that("fixed_d = 1.5 triggers upfront error", {
  expect_error(
    iarimax(dataframe = panel, y_series = "y", x_series = "x",
            id_var = "id", timevar = "time", fixed_d = 1.5),
    regexp = "fixed_d"
  )
})

test_that("fixed_d = 'a' triggers upfront error", {
  expect_error(
    iarimax(dataframe = panel, y_series = "y", x_series = "x",
            id_var = "id", timevar = "time", fixed_d = "a"),
    regexp = "fixed_d"
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
