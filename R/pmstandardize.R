#' Person-level standardize (within-person z-score) time series.
#'
#' Computes within-person z-scores for each variable in `cols`: for every
#' subject, subtracts that subject's mean and divides by that subject's SD.
#'
#' @export
#'
#' @param df A dataframe. Any existing grouping is removed before processing.
#' @param cols A non-empty character vector of column names to standardize.
#' @param id_var A string naming the ID variable (one variable only).
#' @param verbose If `TRUE`, prints a description of the transformation rules
#'   (default `FALSE`).
#' @param append If `TRUE` (default), returns the original dataframe with the
#'   new standardized (`_psd`) columns appended. If `FALSE`, returns only the ID column plus
#'   the new `_psd` columns.
#'
#' @details
#' For each subject, and each column:
#' - All values `NA`: returns `NA`.
#' - Fewer than 2 non-`NA` values: returns `0` (SD undefined; treated as
#'   zero deviation).
#' - Within-person SD ≈ 0 (constant series): returns `0`.
#' - Otherwise: returns `(x - mean(x)) / sd(x)`, computed with `na.rm = TRUE`.
#'
#' Zero-variance subjects are retained as all-zero rows rather than dropped.
#' If you plan to pass output to [iarimax()], those subjects will be filtered
#' out by its `minvar` threshold.
#'
#' @return A dataframe with new `<col>_psd` columns containing the
#'   within-person z-scores.
#'
#' @seealso [i_screener()] for upstream data quality screening,
#'   [i_detrender()] for the next pipeline step, [iarimax()] for per-subject
#'   modeling.
#'
#' @examples
#' local({
#' # Build a small panel: 3 subjects, 10 observations each
#' set.seed(1)
#' panel <- do.call(rbind, lapply(1:3, function(id) {
#'   data.frame(
#'     id   = as.character(id),
#'     time = seq_len(10),
#'     x    = rnorm(10, mean = id * 2),  # different person-means
#'     y    = rnorm(10),
#'     stringsAsFactors = FALSE
#'   )
#' }))
#'
#' # Standardize x and y within each person (append = TRUE, the default)
#' result <- pmstandardize(panel, cols = c("x", "y"), id_var = "id")
#' head(result)  # original columns + x_psd + y_psd
#'
#' # Return only the ID and standardized columns
#' result_slim <- pmstandardize(panel, cols = "x", id_var = "id", append = FALSE)
#' head(result_slim)
#' })

pmstandardize <- function(df, cols, id_var, verbose = FALSE, append = TRUE) {

  # Guard: cols must be non-empty
  if (length(cols) == 0) {
    stop("'cols' must contain at least one column name.")
  }

  # Check if id_var is character and length = 1.
  if (!is.character(id_var) || length(id_var) != 1) {
    stop("'id_var' must be a single character string.")
  }

  # Check if the provided variables are in the dataframe
  required_vars <- c(cols, id_var)

  if (!all(required_vars %in% colnames(df))) {
    missing_vars <- required_vars[!required_vars %in% colnames(df)]
    stop(paste("Cannot find required variables. Check if you spelled the following variables correctly:", paste(missing_vars, collapse = ", ")))
  }

  # Guard: cols must be numeric columns.
  non_numeric <- cols[!vapply(df[cols], is.numeric, logical(1))]
  if (length(non_numeric) > 0) {
    stop(
      "The following columns must be numeric: ",
      paste(non_numeric, collapse = ", "), ". Got class: ",
      paste(vapply(df[non_numeric], function(x) class(x)[1], character(1)), collapse = ", "),
      "."
    )
  }

  # Guard: no Inf/-Inf values in cols (would silently corrupt z-scores).
  has_inf <- cols[vapply(df[cols], function(x) any(is.infinite(x), na.rm = TRUE), logical(1))]
  if (length(has_inf) > 0) {
    stop(
      "Column(s) contain Inf or -Inf values: ",
      paste(has_inf, collapse = ", "),
      ". Remove or replace infinite values before calling pmstandardize()."
    )
  }

  # Remove any pre-existing grouping to avoid silent override
  df <- dplyr::ungroup(df)

  # Provide explanation, conditional to verbose = TRUE.
  if (verbose) {
    message('This function creates within-person z-scores (person-level standardization).')
    message('   If all values for a feature within ID are NA, returns NA.')
    message('   If fewer than 2 non-NA values exist, returns 0 (SD undefined).')
    message('   If values are constant within person, returns 0 (zero variance).')
    message('   Otherwise, returns (x - person_mean) / person_sd.')
  }

  # Apply within-person z-scoring.
  df <- df |>
    dplyr::group_by(!!rlang::sym(id_var)) |>
    dplyr::mutate(dplyr::across(
      dplyr::all_of(cols),
      ~ if (all(is.na(.))) {
          NA_real_
        } else if (sum(!is.na(.)) < 2) {
          ifelse(is.na(.), NA_real_, 0)  # bare 0 would recycle to all positions, overwriting NAs; preserve structure
        } else if (stats::sd(., na.rm = TRUE) < .Machine$double.eps) {
          ifelse(is.na(.), NA_real_, 0)  # same: constant series may still have NA rows; keep them as NA
        } else {
          (. - mean(., na.rm = TRUE)) / stats::sd(., na.rm = TRUE)
        },
      .names = "{.col}_psd"
    )) |>
    dplyr::ungroup()

  # IF append == TRUE, return the original dataframe with new standardized columns appended.
  if (append) {
    return(df)
  }

  # IF append == FALSE, return only id and standardized columns.
  else {
    standardized_cols <- paste0(cols, "_psd")
    return(df |>
      dplyr::select(!!rlang::sym(id_var), dplyr::all_of(standardized_cols)))
  }

}

