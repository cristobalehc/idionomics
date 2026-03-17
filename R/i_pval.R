#' Case by case p-value calculation based on t-distribution.
#'
#' @export
#'
#' @param iarimax_object An iarimax object.
#' @param feature Feature name to calculate p-value. Defaults to iarimax focal predictor. Use your original name, function will automatically append "estimate_"
#' @returns Returns an updated version of results_df within the iarimax_object with p-values for the specific feature.

i_pval <- function(iarimax_object, feature = NULL) {

  #Only works for iarimax object now.
  if (!inherits(iarimax_object, "iarimax_results")) { #will add other classes if needed.
    stop("iarimax_object must be an iarimax_results object.")
  }

  #Defaults to focal predictor if feature is null.
  if (is.null(feature)) {
    feature <- attr(iarimax_object, "focal_predictor")
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
  #Create a mask for valid cases not na and > 0.
  valid_df <- !is.na(df_vec) & df_vec > 0
  #Assign p values only to TRUEs in the mask.
  pval[valid_df] <- 2 * stats::pt(-abs(t_stat[valid_df]), df_vec[valid_df])
  #Assign the pval column to the iarimax_object.
  iarimax_object$results_df[[pval_col_name]] <- pval


  #Return the modified object
  return(invisible(iarimax_object))
}
