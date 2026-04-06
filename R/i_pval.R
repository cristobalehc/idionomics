#' Case by case p-value calculation based on t-distribution.
#'
#' @export
#'
#' @param iarimax_object An iarimax object.
#' @param feature Feature name to calculate p-value. Defaults to iarimax focal predictor. Use your original name, function will automatically append "estimate_"
#' @return Returns an updated version of results_df within the iarimax_object with p-values for the specific feature.
#'   p-values use ML-based degrees of freedom (n_valid - n_params) and will
#'   differ from \code{lm()} even for an ARIMA(0,0,0) model because
#'   \code{stats::arima} estimates residual variance via maximum likelihood
#'   (dividing by n) rather than the OLS unbiased estimator (dividing by
#'   n - k), which affects the standard errors used to compute the t-statistic.
#'   For higher-order models, differencing further reduces n_valid and AR/MA
#'   parameters increase n_params, widening the gap. If both estimate and SE are exactly zero,
#'   the p-value is set to \code{NA} (not \code{NaN}).
#'
#' @examples
#' \donttest{
#' # Fit I-ARIMAX on a small synthetic panel
#' set.seed(42)
#' panel <- do.call(rbind, lapply(1:4, function(id) {
#'   x <- rnorm(30)
#'   data.frame(id = as.character(id), time = seq_len(30),
#'              x = x, y = 0.5 * x + rnorm(30),
#'              stringsAsFactors = FALSE)
#' }))
#'
#' result <- iarimax(panel, y_series = "y", x_series = "x",
#'                   id_var = "id", timevar = "time",
#'                   min_n_subject = 20)
#'
#' # Add per-subject p-values for the focal predictor x
#' result_pval <- i_pval(result)
#' result_pval$results_df[, c("id", "estimate_x", "pval_x")]
#' }

i_pval <- function(iarimax_object, feature = NULL) {

  #Only works for iarimax object.
  if (!inherits(iarimax_object, "iarimax_results")) { #will add other classes if needed.
    stop("iarimax_object must be an iarimax_results object.")
  }

  #Defaults to focal predictor if feature is null.
  if (is.null(feature)) {
    feature <- attr(iarimax_object, "focal_predictor")
  }

  if (!is.character(feature) || length(feature) != 1) {
    stop("'feature' must be a single character string.")
  }

  #Construct the column names based on the 'feature' argument
  feature_name <- paste0("estimate_", feature)
  std_feature_name <- paste0("std.error_", feature)
  pval_col_name <- paste0("pval_",feature)

  #Check if the necessary columns exist
  if (!feature_name %in% names(iarimax_object$results_df)) {
    stop(paste0("Coefficient column '", feature, "' not found in the results data frame."))
  }
  if (!std_feature_name %in% names(iarimax_object$results_df)) {
    stop(paste0("Standard error column '", std_feature_name, "' not found."))
  }

  #calculate df and t statistic.
  df_vec <- iarimax_object$results_df$n_valid - iarimax_object$results_df$n_params
  t_stat <- iarimax_object$results_df[[feature_name]] / iarimax_object$results_df[[std_feature_name]]

  #Check cases where there are negative df's and raise a warning.
  if (any(df_vec <= 0, na.rm = TRUE)) {
    bad_ids <- iarimax_object$results_df[[attr(iarimax_object, "id_var")]][!is.na(df_vec) & df_vec <= 0]
    warning("i_pval: degrees of freedom <= 0 for subject(s): ",
            paste(bad_ids, collapse = ", "),
            ". p-value set to NA. Consider reducing model complexity or increasing min_n_subject.")
  }


  #Create a NA original vector.
  pval <- rep(NA_real_, length(df_vec))
  # Mask: df must be positive AND t_stat must be finite (guards against 0/0 = NaN
  # when both estimate and SE are exactly zero, which would otherwise produce NaN
  # instead of NA).
  valid_df <- !is.na(df_vec) & df_vec > 0 & is.finite(t_stat)
  #Assign p values only to TRUEs in the mask.
  pval[valid_df] <- 2 * stats::pt(-abs(t_stat[valid_df]), df_vec[valid_df])
  #Assign the pval column to the iarimax_object.
  iarimax_object$results_df[[pval_col_name]] <- pval


  #Return the modified object
  return(invisible(iarimax_object))
}
