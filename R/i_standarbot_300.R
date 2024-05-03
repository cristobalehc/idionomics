#' @export
#' @importFrom dplyr %>%

i_standarbot_300 <- function (df, cols, idvar, explanation = TRUE) {

  # df = a dataframe.
  # cols = a vector of strings, with column names.
  # idvar = id variable, just one.

  # This function will return a dataframe with new person standardized columns based on id.
  # If all values for a feature within ID are NA, will return NA.
  # If there is less than two values, it will return 0 (zero variance).
  # If there is no variance, will return 0.
  # Else with return the deviation from the mean.

  if (explanation == TRUE) {
    cat('This function will create a dataframe with person-mean standardized variables appended at the end',"\n")
    Sys.sleep(0.8)
    cat('.')
    Sys.sleep(0.4)
    cat('.')
    Sys.sleep(0.4)
    cat('.',"\n")
    Sys.sleep(0.8)
    cat('If all values for a feature within ID are NA, will return NA',"\n")
    Sys.sleep(0.4)
    cat('If there is less than two values, will return 0 (zero variance)',"\n")
    Sys.sleep(0.4)
    cat('If values are constant, will return 0 (zero variance)',"\n")
    Sys.sleep(0.4)
    cat('If there is enough variance, will return deviations from the mean',"\n")
  }
  #Apply standardization.
  df <- df %>%
    dplyr::group_by(!!rlang::sym(idvar)) %>% #We group by id.
    dplyr::mutate(dplyr::across(dplyr::all_of(cols),  #Mutate across all columns.
                  ~ if (all(is.na(.))) NA_real_ else if (sum(!is.na(.)) <2) 0 else if (sd(., na.rm=TRUE) == 0) 0 else (.-mean(., na.rm=TRUE)) / sd(., na.rm=TRUE), #Apply standardization, catch NA, one value,  and zero variance.
                  .names = "{.col}_PSD")) %>% # Add sufix.
    dplyr::ungroup()

  return(df)

}

