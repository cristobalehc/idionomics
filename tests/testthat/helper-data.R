# Shared helpers loaded automatically by testthat before any test file.

# Build a small reproducible panel dataset.
# Uses arima.sim() with a mild AR(1) so auto.arima reliably selects a
# low-order model, keeping test runtime manageable.
make_panel <- function(n_subjects = 4, n_obs = 25, seed = 42) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(n_subjects), function(id) {
    x <- rnorm(n_obs)
    y <- as.numeric(arima.sim(list(ar = 0.4), n = n_obs)) + 0.5 * x
    data.frame(
      id   = as.character(id),
      time = seq_len(n_obs),
      x    = x,
      y    = y,
      stringsAsFactors = FALSE
    )
  }))
}
