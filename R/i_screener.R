#' Screen subjects for data quality before entering the idionomic pipeline.
#'
#' Applies per-subject data quality filters on raw (unstandardized) time series designed to be implemented
#' before [pmstandardize()] or [iarimax()]. After [pmstandardize()], all
#' non-constant series have within-person variance = 1 by
#' construction, making [iarimax()]'s `minvar` filter ineffective. Running
#' `i_screener()` on raw data prevents low-quality subjects from entering the
#' pipeline.
#'
#' @export
#'
#' @param df A dataframe. Any existing grouping is removed before processing.
#' @param cols A non-empty character vector of column names to screen.
#' @param id_var A string naming the ID variable (one variable only).
#' @param min_n_subject Integer. Subjects or columns with fewer than `min_n_subject`
#'   non-`NA` observations in a given column fail this criterion. Defaults to
#'   `20`, matching [iarimax()]'s default threshold.
#' @param min_sd Numeric or `NULL`. If provided, subjects or columns whose within-person SD
#'   for a given column is below `min_sd` (in raw units) fail this criterion.
#'   Must be a finite positive number. `NULL` (default) skips this check. Use
#'   this before [pmstandardize()] to exclude floor/ceiling responders and
#'   near-constant series.
#' @param max_mode_pct Numeric in `(0, 1]` or `NULL`. If provided, subjects for
#'   whom more than `max_mode_pct` of their non-`NA` responses fall on the same
#'   value fail this criterion. Computed as
#'   `max(table(x)) / sum(!is.na(x))`. Useful for
#'   Likert/ordinal EMA data to detect "stuck" or floor/ceiling responders.
#'   `NULL` (default) skips this check.
#' @param filter_type One of `"joint"` (default) or `"per_column"`. Controls
#'   how criteria are applied across multiple columns:
#'   - `"joint"`: a subject is excluded or flagged if they fail *any* criterion
#'     on *any* variable. Mirrors [iarimax()]'s AND filter. Recommended when
#'     all `cols` will be analyzed together in the same model.
#'   - `"per_column"`: each column is evaluated independently. A subject can
#'     pass for one variable and fail for another. Useful for exploratory quality
#'     inspection across a broad set of variables. The pass_overall column in reports always uses a joint criterion.
#' @param mode One of `"filter"` (default), `"flag"`, or `"report"`. Controls
#'   the type of output returned:
#'   - `"filter"`: returns the dataframe with failing subjects' rows removed
#'     (`"joint"`) or their failing column values set to `NA` (`"per_column"`).
#'   - `"flag"`: returns the original dataframe with a logical `pass_overall`
#'     column added (`"joint"`) or one `<col>_pass` column per variable
#'     (`"per_column"`).
#'   - `"report"`: returns a per-subject summary dataframe with quality metrics
#'     and pass/fail indicators for each column and overall (pass_overall).
#' @param verbose If `TRUE`, prints the criteria applied and subject counts
#'   before and after screening. Defaults to `FALSE`.
#'
#' @details
#' For each subject x column combination, three quality metrics are computed:
#' - `n_valid`: number of non-`NA` observations.
#' - `sd`: within-person standard deviation (with `na.rm = TRUE`); `NA` when
#'   fewer than 2 non-`NA` values exist.
#' - `mode_pct`: proportion of non-`NA` responses equal to the modal value,
#'   computed as `max(table(x)) / length(x[!is.na(x)])`. Uses `table()`
#'   internally, so values are compared as characters; best suited to
#'   integer or Likert-scale data where floating-point rounding is not a concern.
#'
#' Criteria are applied in order: `min_n_subject`, then `min_sd`, then `max_mode_pct`.
#' A subject-column fails as soon as any active criterion is not met.
#'
#' When `filter_type = "per_column"` and `mode = "filter"`, failing subjects
#' are not removed entirely — their values in the failing column are set to
#' `NA`. Subject-column combinations consisting solely of NA values will fail downstream function's `min_n_subject` filters.
#'
#' Note that `i_screener()` evaluates each column's non-`NA` count independently
#' (`min_n`), whereas [iarimax()] filters on *pairwise-complete* observations
#' across all series jointly. A subject that passes `i_screener()`'s `min_n`
#' threshold may still be excluded by [iarimax()] if the non-`NA` rows in the
#' outcome and predictor series do not sufficiently overlap.
#'
#' When `mode = "report"`, `filter_type` does not affect the output shape — the
#' report always contains one row per subject with quality metrics for all
#' columns. `pass_overall` always reflects the joint AND across all columns,
#' matching `filter_type = "joint"` semantics.
#'
#' The `minvar` filters in [iarimax()] and [i_detrender()] serve a different,
#' complementary role: they are technical last-resort guards that protect each
#' function when called independently. In particular, [i_detrender()]'s
#' post-detrend variance check catches series that pass `i_screener()` on raw
#' data but produce near-zero residuals after removing a near-perfect linear
#' trend — a case `i_screener()` cannot anticipate.
#'
#' @return Depends on `mode`:
#' - `"filter"`: a dataframe with the same columns as input but possibly fewer
#'   rows (`"joint"`) or `NA` values in screened columns (`"per_column"`).
#' - `"flag"`: the original dataframe with `pass_overall` (`"joint"`) or
#'   `<col>_pass` columns (`"per_column"`) appended.
#' - `"report"`: a per-subject summary dataframe with one row per subject and
#'   columns grouped by metric: all `<col>_n_valid`, then all `<col>_sd`, then
#'   all `<col>_mode_pct`, then all `<col>_pass`, and finally `pass_overall`.
#'
#' @examples
#' local({
#' set.seed(42)
#' panel <- do.call(rbind, lapply(1:4, function(id) {
#'   data.frame(
#'     id   = as.character(id),
#'     time = seq_len(25),
#'     x    = if (id == 2) rep(3L, 25) else sample(1:7, 25, replace = TRUE),
#'     y    = sample(1:7, 25, replace = TRUE),
#'     stringsAsFactors = FALSE
#'   )
#' }))
#'
#' # Remove subjects failing default min_n_subject = 20 or a minimum SD criterion
#' result <- i_screener(panel, cols = c("x", "y"), id_var = "id", min_sd = 0.5)
#'
#' # Inspect which subjects would be removed without committing
#' flagged <- i_screener(panel, cols = c("x", "y"), id_var = "id",
#'                     min_sd = 0.5, mode = "flag")
#' table(flagged$pass_overall)
#'
#' # Retrieve a per-subject quality summary
#' report <- i_screener(panel, cols = c("x", "y"), id_var = "id",
#'                    min_sd = 0.5, mode = "report")
#' print(report)
#' })

i_screener <- function(df, cols, id_var,
                     min_n_subject = 20,
                     min_sd        = NULL,
                     max_mode_pct = NULL,
                     filter_type  = "joint",
                     mode         = "filter",
                     verbose      = FALSE) {

  # Guard: cols must be non-empty.
  if (length(cols) == 0) {
    stop("'cols' must contain at least one column name.")
  }

  if (!is.character(id_var) || length(id_var) != 1) {
    stop("'id_var' must be a single character string.")
  }

  # Check if the provided variables are in the dataframe.
  required_vars <- c(cols, id_var)
  if (!all(required_vars %in% colnames(df))) {
    missing_vars <- required_vars[!required_vars %in% colnames(df)]
    stop(paste(
      "Cannot find required variables. Check if you spelled the following variables correctly:",
      paste(missing_vars, collapse = ", ")
    ))
  }

  # Validate threshold parameters.
  if (!is.numeric(min_n_subject) || length(min_n_subject) != 1 ||
      !is.finite(min_n_subject) || min_n_subject < 1) {
    stop("'min_n_subject' must be a finite positive number.")
  }
  if (!is.null(min_sd) && (!is.numeric(min_sd) || length(min_sd) != 1 ||
      !is.finite(min_sd) || min_sd <= 0)) {
    stop("'min_sd' must be a finite positive number.")
  }
  if (!is.null(max_mode_pct) && (!is.numeric(max_mode_pct) || length(max_mode_pct) != 1 ||
      max_mode_pct <= 0 || max_mode_pct > 1)) {
    stop("'max_mode_pct' must be a number strictly greater than 0 and at most 1.")
  }
  if (!filter_type %in% c("joint", "per_column")) {
    stop("'filter_type' must be 'joint' or 'per_column'.")
  }
  if (!mode %in% c("filter", "flag", "report")) {
    stop("'mode' must be 'filter', 'flag', or 'report'.")
  }

  # Check for output column name collisions that would corrupt output silently.
  if (mode == "flag" && filter_type == "joint" && "pass_overall" %in% names(df)) {
    stop(
      "Input dataframe already contains a column named 'pass_overall'. ",
      "Please rename it before calling i_screener()."
    )
  }
  if (filter_type == "per_column") {
    conflicts <- intersect(names(df), paste0(cols, "_pass"))
    if (length(conflicts) > 0) {
      stop(
        "Input dataframe already contains column(s) that i_screener uses internally: ",
        paste(conflicts, collapse = ", "),
        ". Please rename them before calling i_screener()."
      )
    }
  }

  # Remove any pre-existing grouping.
  df <- dplyr::ungroup(df)

  #Create symbolic variable.
  id_var_sym <- rlang::sym(id_var)

  #Calculate the original n_subjects.
  n_subjects_original <- length(unique(df[[id_var]]))

  # Provide explanation, conditional to verbose = TRUE.
  if (verbose) {
    message("i_screener applies per-subject data quality filters.")
    message("   min_n_subject: subject-column combinations need >= ", min_n_subject, " non-NA observations per variable.")
    if (!is.null(min_sd)) {
      message("   min_sd       : subject-column combinations need within-person SD >= ", min_sd, " per variable.")
    }
    if (!is.null(max_mode_pct)) {
      message("   max_mode_pct : subject-column combinations need <= ", max_mode_pct * 100, "% of responses on the modal value.")
    }
    message(
      "   filter_type  : '", filter_type, "' - ",
      if (filter_type == "joint") "Subjects excluded if any variable fails any criterion."
      else "Each variable evaluated independently."
    )
  }

  # Compute per-subject, per-column quality metrics.
  metrics <- df |>
    dplyr::group_by(!!id_var_sym) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(cols),
        list(
          n_valid  = ~ sum(!is.na(.x)),
          # sd() requires >= 2 non-NA values to return a finite result.
          sd       = ~ if (sum(!is.na(.x)) >= 2) stats::sd(.x, na.rm = TRUE)
                       else NA_real_,
          mode_pct = ~ {
            x_obs <- .x[!is.na(.x)]
            if (length(x_obs) == 0) NA_real_
            else as.numeric(max(table(x_obs))) / length(x_obs)
          }
        ),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    )

  # Apply each criterion and record pass/fail per subject per column.
  for (col in cols) {
    n_col    <- paste0(col, "_n_valid")
    sd_col   <- paste0(col, "_sd")
    mode_col <- paste0(col, "_mode_pct")
    pass_col <- paste0(col, "_pass")

    #Per subject pass: min_n_subject.
    col_pass <- metrics[[n_col]] >= min_n_subject

    if (!is.null(min_sd)) {
      col_pass <- col_pass & !is.na(metrics[[sd_col]]) & (metrics[[sd_col]] >= min_sd)
    }
    if (!is.null(max_mode_pct)) {
      col_pass <- col_pass & !is.na(metrics[[mode_col]]) & (metrics[[mode_col]] <= max_mode_pct)
    }

    metrics[[pass_col]] <- col_pass
  }

  # Compute overall pass: subject passes only if every column passes every criterion.
  pass_cols <- paste0(cols, "_pass")
  metrics <- metrics |>
    dplyr::mutate(
      pass_overall = dplyr::if_all(dplyr::all_of(pass_cols))
    )

  if (verbose) {
    n_pass <- sum(metrics$pass_overall)
    n_fail <- n_subjects_original - n_pass
    if (length(cols) == 1) {
      fail_label <- "Subjects failing at least one criterion"
      pass_label <- "Subjects passing all criteria"
    } else {
      fail_label <- "Subjects failing at least one criterion across variables"
      pass_label <- "Subjects passing all criteria on all variables"
    }
    message("   Subjects before screening : ", n_subjects_original)
    message("   ", fail_label, " : ", n_fail)
    message("   ", pass_label, " : ", n_pass)
  }

  # IF mode == "report", return a per-subject quality summary table.
  if (mode == "report") {
    out_cols <- c( #Create vector with names, to reorder output by metric.
      id_var,
      paste0(cols, "_n_valid"),
      paste0(cols, "_sd"),
      paste0(cols, "_mode_pct"),
      paste0(cols, "_pass"),
      "pass_overall"
    )
    return(metrics |> dplyr::select(dplyr::all_of(out_cols)))
  }

  # IF mode == "flag", append pass/fail columns to the original dataframe.
  if (mode == "flag") {
    if (filter_type == "joint") {
      pass_info <- metrics[, id_var, drop = FALSE]
      pass_info$pass_overall <- metrics$pass_overall
      return(dplyr::left_join(df, pass_info, by = id_var))
    } #if filter_type = "per_column"
    else {
      pass_info <- metrics |>
        dplyr::select(!!id_var_sym, dplyr::all_of(pass_cols))
      return(dplyr::left_join(df, pass_info, by = id_var))
    }
  }
  #Else:  mode = filter.
  # IF mode == "filter" + filter_type == "joint", remove failing subjects entirely.
  # IF mode == "filter" + filter_type == "per_column", set failing column values to NA.
  if (filter_type == "joint") {
    keep_ids <- metrics[[id_var]][metrics$pass_overall]
    return(df |> dplyr::filter(!!id_var_sym %in% keep_ids))
  }
  else {
    pass_info <- metrics |>
      dplyr::select(!!id_var_sym, dplyr::all_of(pass_cols))
    df <- dplyr::left_join(df, pass_info, by = id_var)
    #Set failing column values to NA.
    for (col in cols) {
      pass_col <- paste0(col, "_pass")
      df[[col]][!df[[pass_col]]] <- NA #For each element in col, set those not pass (false) as NA.
      df[[pass_col]] <- NULL
    }
    return(df)
  }

}
