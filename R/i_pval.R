#' Case by case p-value calculation based on t-distribution.
#'
#' @export
#'
#' @param iarimax_object An iarimax object.
#' @param feature Feature name to calculate p-value. Defaults to "xreg".
#' @returns Returns an updated version of results_df within the iarimax_object with p-values for the specific feature.

i_pval <- function(iarimax_object, feature = "xreg") {

  #Construct the column names based on the 'feature' argument
  std_feature_name <- paste0("stderr_", feature)
  pval_col_name <- paste0("pval_",feature)

  #Check if the necessary columns exist
  if (!feature %in% names(iarimax_object$results_df)) {
    stop(paste0("Error: Coefficient column '", feature, "' not found in the results data frame."))
  }
  if (!std_feature_name %in% names(iarimax_object$results_df)) {
    stop(paste0("Error: Standard error column '", std_feature_name, "' not found."))
  }

  #Calculate p_value
  iarimax_object$results_df[[pval_col_name]] <- 2*stats::pt(-abs(iarimax_object$results_df[[feature]] /
                                                                   iarimax_object$results_df[[std_feature_name]]),
                                                            iarimax_object$results_df$n_valid - iarimax_object$results_df$n_params)


  #Return the modified object
  return(iarimax_object)
}
