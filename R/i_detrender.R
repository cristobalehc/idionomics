#' Detrend idiographic time series by time variable.
#'
#' @export
#' @importFrom dplyr %>%
#'
#' @param df A dataframe.
#' @param cols A vector of strings with the variable names you would like to detrend.
#' @param idvar Your ID variable. It currently accepts only one.
#' @param timevar The time variable used to detrend yout data.
#' @param explanation Boolean to print explanation of the function (default is TRUE).
#' @param append Boolean to append new standardized columns to your original dataframe or not (default is TRUE)
#'
#' @returns Will return a dataframe with new detrended columns appended at the end.

i_detrender <- function (df, cols, idvar, timevar, explanation = TRUE, append = TRUE) {

  # df = a dataframe.
  # cols = a vector of strings, with column names.
  # idvar = id variable, just one.
  # timevar = a time variable.

  # This function will return a dataframe with detrended timeseries per individual.
  # If all values for a feature within ID are NA, will return NA.
  # If there is less than three values, it will return NA.
  # Else with return the residuals from the linear model with time.

  # Check if the provided variables are in the dataframe
  required_vars <- c(cols, idvar, timevar)

  if (!all(required_vars %in% colnames(df))) {
    missing_vars <- required_vars[!required_vars %in% colnames(df)]
    stop(paste("Cannot find required variables. Check if you spelled the following variables correctly:", paste(missing_vars, collapse = ", ")))
  }


  # Provide explanation, conditional to explanation = TRUE.
  if (explanation == TRUE) {
    cat('This function will return a dataframe with detrended timeseries per individual.',"\n")
    Sys.sleep(0.8)
    cat('.')
    Sys.sleep(0.4)
    cat('.')
    Sys.sleep(0.4)
    cat('.',"\n")
    Sys.sleep(0.8)
    cat('   If all values for a feature within ID are NA, will return NA',"\n")
    Sys.sleep(0.4)
    cat('   If there is less than three values (perfect line, zero residual), will return NA',"\n")
    Sys.sleep(0.4)
    cat('   Else, will return the residuals from the linear model with time.',"\n")
  }

  #Detrend

  df <- df %>%
    dplyr::group_by(!!rlang::sym(idvar)) %>% #We group by id.
    dplyr::mutate(dplyr::across(dplyr::all_of(cols),  #Mutate across all columns.
                                ~ if (all(is.na(.))) NA_real_ else if (sum(!is.na(.)) <3) NA_real_ else stats::residuals(stats::lm(. ~ .data[[timevar]], na.action = stats::na.exclude)), #Detrend
                                .names = "{.col}_DT")) %>% # Add sufix "DT" for "DeTrended"
    dplyr::ungroup()


  # IF append == TRUE, return the original dataframe with new detrended columns appended at the end (right)
  if (append == TRUE) {

    cat("\n")
    cat('   Detrend finished: Your original dataframe with detrended columns appended at the end was returned')
    return(df) }

  # IF append == FALSE, return a smaller dataframe only with id, timevar and detrended columns.
  else {
    #Paste column names with _PSD
    detrended_cols = paste0(cols, "_DT")

    #Filter those columns only.
    df <- df %>%
      dplyr::select(!!rlang::sym(idvar), !!rlang::sym(timevar), dplyr::all_of(detrended_cols))

    cat('   Detrend finished: A dataframe only with your ID, time variable and detrended columns was returned')
    return(df)
  }

}

