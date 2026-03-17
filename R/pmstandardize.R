#' Person-level standardize (within-person z-score) time series.
#'
#' Computes within-person z-scores for each variable in `cols`: for every
#' subject, subtracts that subject's mean and divides by that subject's SD.
#' This is stronger than person-mean centering (mean subtraction only) — it
#' also removes between-person differences in variance, making within-person
#' fluctuations unit-free and comparable across subjects.
#'
#' @export
#'
#' @param df A dataframe. Must not already be grouped; any existing grouping is
#'   removed before processing.
#' @param cols A non-empty character vector of column names to standardize.
#' @param idvar A string naming the ID variable (one variable only).
#' @param verbose If `TRUE`, prints a description of the transformation rules
#'   (default `FALSE`).
#' @param append If `TRUE` (default), returns the original dataframe with the
#'   new `_psd` columns appended. If `FALSE`, returns only the ID column plus
#'   the new `_psd` columns.
#'
#' @details
#' For each subject × column combination:
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
#' @returns A dataframe with new `<col>_psd` columns containing the
#'   within-person z-scores.

pmstandardize <- function(df, cols, idvar, verbose = FALSE, append = TRUE) {

  # Guard: cols must be non-empty
  if (length(cols) == 0) {
    stop("'cols' must contain at least one column name.")
  }

  # Check if the provided variables are in the dataframe
  required_vars <- c(cols, idvar)

  if (!all(required_vars %in% colnames(df))) {
    missing_vars <- required_vars[!required_vars %in% colnames(df)]
    stop(paste("Cannot find required variables. Check if you spelled the following variables correctly:", paste(missing_vars, collapse = ", ")))
  }

  # Remove any pre-existing grouping to avoid silent override
  df <- dplyr::ungroup(df)

  # Provide explanation, conditional to verbose = TRUE.
  if (verbose == TRUE) {
    message('This function creates within-person z-scores (person-level standardization).', "\n")
    message('   If all values for a feature within ID are NA, returns NA.', "\n")
    message('   If fewer than 2 non-NA values exist, returns 0 (SD undefined).', "\n")
    message('   If values are constant within person, returns 0 (zero variance).', "\n")
    message('   Otherwise, returns (x - person_mean) / person_sd.', "\n")
  }

  # Apply within-person z-scoring.
  df <- df |>
    dplyr::group_by(!!rlang::sym(idvar)) |>
    dplyr::mutate(dplyr::across(
      dplyr::all_of(cols),
      ~ if (all(is.na(.))) {
          NA_real_
        } else if (sum(!is.na(.)) < 2) {
          ifelse(is.na(.), NA_real_, 0)  # bare 0 would recycle to all positions, overwriting NAs; preserve structure
        } else if (sd(., na.rm = TRUE) < .Machine$double.eps) {
          ifelse(is.na(.), NA_real_, 0)  # same: constant series may still have NA rows; keep them as NA
        } else {
          (. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE)
        },
      .names = "{.col}_psd"
    )) |>
    dplyr::ungroup()

  # IF append == TRUE, return the original dataframe with new standardized columns appended.
  if (append == TRUE) {
    return(df)
  }

  # IF append == FALSE, return only id and standardized columns.
  else {
    standardized_cols <- paste0(cols, "_psd")
    return(df |>
      dplyr::select(!!rlang::sym(idvar), dplyr::all_of(standardized_cols)))
  }

}
