#' Person mean standardize time series.
#'
#' @export
#' @importFrom dplyr %>%
#'
#' @param df A dataframe.
#' @param cols A vector of strings with the variable names you would like to standardize.
#' @param idvar Your ID variable. It currently accepts only one.
#' @param explanation Boolean to print explanation of the function (default is TRUE).
#' @param append Boolean to append new standardized columns to your original dataframe or not (default is TRUE)
#'
#' @returns Will return a dataframe with new person-mean standardized columns appended at the end.

i_standarbot_300 <- function (df, cols, idvar, explanation = TRUE, append = TRUE) {

  # df = a dataframe.
  # cols = a vector of strings, with column names.
  # idvar = id variable, just one.

  # This function will return a dataframe with new person standardized columns based on id.
  # If all values for a feature within ID are NA, will return NA.
  # If there is less than two values, it will return 0 (zero variance).
  # If there is no variance, will return 0.
  # Else with return the deviation from the mean.

  # Check if the provided variables are in the dataframe
  required_vars <- c(cols, idvar)

  if (!all(required_vars %in% colnames(df))) {
    missing_vars <- required_vars[!required_vars %in% colnames(df)]
    stop(paste("Cannot find required variables. Check if you spelled the following variables correctly:", paste(missing_vars, collapse = ", ")))
  }


  # Provide explanation, conditional to explanation = TRUE.
  if (explanation == TRUE) {
    cat('This function will create a dataframe with person-mean standardized variables appended at the end',"\n")
    Sys.sleep(0.8)
    cat('.')
    Sys.sleep(0.4)
    cat('.')
    Sys.sleep(0.4)
    cat('.',"\n")
    Sys.sleep(0.8)
    cat('   If all values for a feature within ID are NA, will return NA',"\n")
    Sys.sleep(0.4)
    cat('   If there is less than two values, will return 0 (zero variance)',"\n")
    Sys.sleep(0.4)
    cat('   If values are constant, will return 0 (zero variance)',"\n")
    Sys.sleep(0.4)
    cat('   If there is enough variance, will return deviations from the mean',"\n")
  }

  #Apply standardization.
  df <- df %>%
    dplyr::group_by(!!rlang::sym(idvar)) %>% #We group by id.
    dplyr::mutate(dplyr::across(dplyr::all_of(cols),  #Mutate across all columns.
                  ~ if (all(is.na(.))) NA_real_ else if (sum(!is.na(.)) <2) 0 else if (sd(., na.rm=TRUE) == 0) 0 else (.-mean(., na.rm=TRUE)) / sd(., na.rm=TRUE), #Apply standardization, catch NA, one value,  and zero variance.
                  .names = "{.col}_PSD")) %>% # Add sufix.
    dplyr::ungroup() %>%
    dplyr::filter()


  # IF append == TRUE, return the original dataframe with new standardized columns appended at the end (right)
  if (append == TRUE) {

  cat("\n")
    cat('   Standardization finished: Your original dataframe with standardized columns appended at the end was returned')
  return(df) }

  # IF append == FALSE, return a smaller dataframe only with id and standardized columns.
  else {
    #Paste column names with _PSD
    standardized_cols = paste0(cols, "_PSD")

    #Filter those columns only.
    df <- df %>%
      dplyr::select(!!rlang::sym(idvar), dplyr::all_of(standardized_cols))

    cat('   Standardization finished: A dataframe only with your ID and standardized columns was returned')
    return(df)
  }

}

