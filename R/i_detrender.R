#' Linearly detrend idiographic time series by the time variable.
#'
#' @export
#'
#' @param df A dataframe. Any existing grouping is removed before processing.
#' @param cols A non-empty character vector of column names to detrend.
#' @param id_var A string naming the ID variable (one variable only).
#' @param timevar A string naming the time variable used for linear detrending. Must
#'   be numeric and complete — no missing values allowed.
#' @param min_n_subject Integer. Subjects with fewer than `min_n_subject`
#'   non-`NA` observations in a given column receive `NA` in the detrended
#'   output. Mirrors the threshold used in [iarimax()]. Defaults to 20.
#' @param minvar Numeric. Last-resort guard against near-zero variance: subjects
#'   whose pre-detrend or post-detrend variance for a given column falls below
#'   `minvar` receive `NA` in the detrended output. Defaults to 0.01. This guard
#'   protects against constant or near-constant series that would produce
#'   meaningless residuals. For substantive data quality screening on raw data,
#'   use [i_screener()] before entering the pipeline.
#' @param verbose If `TRUE`, prints a description of the detrending rules.
#'   Defaults to `FALSE`.
#' @param append If `TRUE` (default), returns the original dataframe with new
#'   `<col>_dt` columns appended. If `FALSE`, returns only the ID column, the
#'   time variable, and the new `_dt` columns.
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
#' a subject can receive `NA` in one `_dt` column while still producing valid
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
#' @return A dataframe with new `<col>_dt` columns containing the
#'   within-person linear-detrended residuals.
#'
#' @examples
#' local({
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
#' # Detrend x within each person (returns original df + x_dt)
#' result <- i_detrender(panel, cols = "x", id_var = "id", timevar = "time")
#' head(result)
#'
#' # Return only id, time, and detrended columns
#' result_slim <- i_detrender(panel, cols = "x", id_var = "id",
#'                            timevar = "time", append = FALSE)
#' head(result_slim)
#' })

i_detrender <- function(df, cols, id_var, timevar,
                        min_n_subject = 20, minvar = 0.01,
                        verbose = FALSE, append = TRUE) {

  # Guard: cols must be non-empty
  if (length(cols) == 0) {
    stop("'cols' must contain at least one column name.")
  }

  #id var as single character string.
  if (!is.character(id_var) || length(id_var) != 1) {
    stop("'id_var' must be a single character string.")
  }

  #timevar single character string.
  if (!is.character(timevar) || length(timevar) != 1) {
    stop("'timevar' must be a single character string.")
  }

  #min_n_subject as finite positive number.
  if (!is.numeric(min_n_subject) || length(min_n_subject) != 1 ||
      !is.finite(min_n_subject) || min_n_subject < 1) {
    stop("'min_n_subject' must be a finite positive number.")
  }

  #minvar finite non-negative number.
  if (!is.numeric(minvar) || length(minvar) != 1 ||
      !is.finite(minvar) || minvar < 0) {
    stop("'minvar' must be a finite non-negative number.")
  }

  # Check if the provided variables are in the dataframe
  required_vars <- c(cols, id_var, timevar)

  if (!all(required_vars %in% colnames(df))) {
    missing_vars <- required_vars[!required_vars %in% colnames(df)]
    stop(paste(
      "Cannot find required variables. Check if you spelled the following variables correctly:",
      paste(missing_vars, collapse = ", ")
    ))
  }

  # timevar must be numeric for a meaningful linear trend
  if (!is.numeric(df[[timevar]])) {
    stop(
      "Column '", timevar, "' must be numeric for linear detrending. Got class: ",
      class(df[[timevar]])[1], ". Convert dates or other formats to a reasonable numeric value before calling i_detrender."
    )
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

  #Create symbolic variables.
  id_var_sym   <- rlang::sym(id_var)
  timevar_sym <- rlang::sym(timevar)

  #Run linear detrending.
  df <- df |>
    dplyr::group_by(!!id_var_sym) |> #Group by id variable.
    dplyr::mutate(dplyr::across( #loop through all columns.
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
      .names = "{.col}_dt"
    )) |>
    dplyr::ungroup()

  if (append) {
    return(df)
  } else {
    dt_cols <- paste0(cols, "_dt")
    return(df |>
      dplyr::select(!!id_var_sym, !!timevar_sym, dplyr::all_of(dt_cols)))
  }
}
