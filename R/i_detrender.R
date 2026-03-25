#' Detrend idiographic time series by the time variable.
#'
#' @export
#'
#' @param df A dataframe. Any existing grouping is removed before processing.
#' @param cols A non-empty character vector of column names to detrend.
#' @param idvar A string naming the ID variable (one variable only).
#' @param timevar A string naming the time variable used for detrending. Must
#'   be numeric (or coercible to numeric) and complete — no missing values
#'   allowed.
#' @param min_n_subject Integer. Subjects with fewer than `min_n_subject`
#'   non-`NA` observations in a given column receive `NA` in the detrended
#'   output. Mirrors the threshold used in [iarimax()]. Defaults to 20.
#' @param minvar Numeric. Subjects whose pre-detrend or post-detrend variance
#'   for a given column falls below `minvar` receive `NA` in the detrended
#'   output. Mirrors the threshold used in [iarimax()]. Defaults to 0.01.
#' @param verbose If `TRUE`, prints a description of the detrending rules.
#'   Defaults to `FALSE`.
#' @param append If `TRUE` (default), returns the original dataframe with new
#'   `<col>_DT` columns appended. If `FALSE`, returns only the ID column, the
#'   time variable, and the new `_DT` columns.
#'
#' @details
#' For each subject x column combination, the following rules are applied in
#' order:
#'
#' - All values `NA`: returns `NA`.
#' - Fewer than `min_n_subject` non-`NA` values: returns `NA`.
#' - Pre-detrend variance < `minvar`: returns `NA`.
#' - Otherwise: fits `lm(col ~ timevar, na.action = na.exclude)` and takes the
#'   residuals. If the post-detrend residual variance is also < `minvar`,
#'   returns `NA`.
#'
#' The thresholds `min_n_subject` and `minvar` share the same defaults as
#' [iarimax()], but filtering here is applied **per column independently**:
#' a subject can receive `NA` in one `_DT` column while still producing valid
#' residuals in another. This differs from [iarimax()], which applies a joint
#' AND filter — a subject is excluded if *any* series fails the thresholds.
#' The per-column design is intentional: it avoids unnecessary data loss and
#' lets you detrend variables in different combinations before deciding which
#' to pass to [iarimax()]. The thresholds still serve the same protective
#' purpose — preventing data-poor or low-variance subjects from being detrended,
#' which would amplify numerical noise to apparent unit variance after
#' person-mean standardization and allow them to slip through [iarimax()]'s
#' filters.
#'
#' @return A dataframe with new `<col>_DT` columns containing the
#'   within-person linear-detrended residuals.
#'
#' @examples
#' # Build a panel with a linear time trend embedded in x
#' set.seed(2)
#' panel <- do.call(rbind, lapply(1:3, function(id) {
#'   data.frame(
#'     id   = as.character(id),
#'     time = seq_len(25),
#'     x    = seq_len(25) + rnorm(25),  # linear trend + noise
#'     stringsAsFactors = FALSE
#'   )
#' }))
#'
#' # Detrend x within each person (returns original df + x_DT)
#' result <- i_detrender(panel, cols = "x", idvar = "id", timevar = "time")
#' head(result)
#'
#' # Return only id, time, and detrended columns
#' result_slim <- i_detrender(panel, cols = "x", idvar = "id",
#'                            timevar = "time", append = FALSE)
#' head(result_slim)

i_detrender <- function(df, cols, idvar, timevar,
                        min_n_subject = 20, minvar = 0.01,
                        verbose = FALSE, append = TRUE) {

  # Guard: cols must be non-empty
  if (length(cols) == 0) {
    stop("'cols' must contain at least one column name.")
  }

  # Check if the provided variables are in the dataframe
  required_vars <- c(cols, idvar, timevar)

  if (!all(required_vars %in% colnames(df))) {
    missing_vars <- required_vars[!required_vars %in% colnames(df)]
    stop(paste(
      "Cannot find required variables. Check if you spelled the following variables correctly:",
      paste(missing_vars, collapse = ", ")
    ))
  }

  # Check for NA in timevar: required for correct temporal ordering in lm
  n_missing_timevar <- sum(is.na(df[[timevar]]))
  if (n_missing_timevar > 0) {
    stop(
      n_missing_timevar, " row(s) have missing values in the time variable '", timevar, "'. ",
      "i_detrender requires a non-missing timevar for every observation. ",
      "Please review, remove or impute these rows before running i_detrender."
    )
  }

  # Remove any pre-existing grouping
  df <- dplyr::ungroup(df)

  if (verbose) {
    message("i_detrender removes the within-person linear time trend per variable.")
    message("   If all values for a column within ID are NA: returns NA.")
    message("   If fewer than ", min_n_subject, " non-NA observations: returns NA.")
    message("   If pre-detrend or post-detrend variance < ", minvar, ": returns NA.")
    message("   Otherwise: returns residuals from lm(col ~ timevar) per person.")
  }

  idvar_sym   <- rlang::sym(idvar)
  timevar_sym <- rlang::sym(timevar)

  df <- df |>
    dplyr::group_by(!!idvar_sym) |>
    dplyr::mutate(dplyr::across(
      dplyr::all_of(cols),
      ~ {
        n_obs   <- sum(!is.na(.))
        pre_var <- stats::var(., na.rm = TRUE)

        if (all(is.na(.))) {
          NA_real_
        } else if (n_obs < min_n_subject || is.na(pre_var) || pre_var < minvar) {
          NA_real_
        } else {
          resid_vals <- stats::residuals(
            stats::lm(. ~ .data[[timevar]], na.action = stats::na.exclude)
          )
          post_var <- stats::var(resid_vals, na.rm = TRUE)
          if (is.na(post_var) || post_var < minvar) {
            NA_real_
          } else {
            resid_vals
          }
        }
      },
      .names = "{.col}_DT"
    )) |>
    dplyr::ungroup()

  if (append) {
    return(df)
  } else {
    dt_cols <- paste0(cols, "_DT")
    return(df |>
      dplyr::select(!!idvar_sym, !!timevar_sym, dplyr::all_of(dt_cols)))
  }
}
