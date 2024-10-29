#' Easy copyMeans imputation for time-series.
#'
#' @export
#' @importFrom dplyr %>%
#' @importFrom rlang `:=`
#'
#' @param data A string containing your dataframe name.
#' @param id_col  A string containing your ID column.
#' @param time_col  A string containing your time column.
#' @param cols_of_interest A vector of strings containing the columns you wish to impute.
#'
#' @returns It returns a list with a new dataframe with imputed values, longdata objects with and without imputation, your filtered data, wide data used to create the longdata object, and a list with indexes to create the longdata object.


# Copy means imputation function:  ~~ IMPUTATRON ~~
imputatron_2000 <- function(data, id_col, time_col, cols_of_interest) {

  #Documentation: This function will take a dataframe in long format, with an id column, a time column and columns of interest
  #               and will produce a copyMeans imputation for every time series for each subject.

  #Assumptions: It assumes that you have already filtered your data (e.g.: selected cases with less than n NA values)
  #             It also assumes that the data is in the long format where each observation is in a row, and variables are in columns.

  cat('~~ Beep Boop Beep Boop ~~',"\n")
  Sys.sleep(0.4)
  cat('Imputatron will try to impute your dataset',"\n")
  Sys.sleep(0.4)
  cat('.')
  Sys.sleep(0.4)
  cat('.')
  Sys.sleep(0.4)
  cat('.',"\n")
  Sys.sleep(0.4)


  # Filter and arrange data
  filtered_data <- data %>%
    dplyr::select(dplyr::all_of(id_col), dplyr::all_of(time_col), dplyr::all_of(cols_of_interest)) %>%
    dplyr::group_by(!!rlang::sym(id_col)) %>%
    dplyr::arrange(!!rlang::sym(id_col), !!rlang::sym(time_col))
  cat('Imputatron filtered your data',"\n")
  Sys.sleep(0.8)
  # Convert to wide format
  wide_data <- filtered_data %>%
    tidyr::pivot_wider(names_from = !!rlang::sym(time_col), values_from = dplyr::all_of(cols_of_interest))
  cat('Imputatron transformed data to wide format',"\n")
  Sys.sleep(0.8)

  # Generate column ranges (it assumes an ID col at position 1 and wide format)
  column_ranges <- lapply(seq_along(cols_of_interest), function(i) {
    start_index <- (i - 1) * max(filtered_data[[rlang::sym(time_col)]]) + 2  # Assuming the first column might be an ID or other non-variable column
    end_index <- i * max(filtered_data[[rlang::sym(time_col)]]) + 1 #It assumes that every variable has the same number of measdurements (for max filtered_data$time_col)
    return(start_index:end_index)
  })
  cat('Imputatron extracted start and end indexes for each timeseries',"\n")
  Sys.sleep(0.8)
  # Create the timeInData list
  timeInData_list <- stats::setNames(column_ranges, cols_of_interest)
  cat('Imputatron dynamically created a list of variables + indexes',"\n")
  Sys.sleep(0.8)
  cat('.')
  Sys.sleep(0.4)
  cat('.')
  Sys.sleep(0.4)
  cat('.',"\n")
  Sys.sleep(0.4)

  cat('###### Imputatron completed data preprocessing ####',"\n")
  Sys.sleep(0.8)
  # Create longitudinal data object
  long_data_obj <- longitudinalData::longData3d(
    wide_data,
    idAll = wide_data[[id_col]],
    timeInData = timeInData_list)
  cat('Imputatron created longData3d object',"\n")
  Sys.sleep(0.8)
  # Perform imputation
  imputed_data <- longitudinalData::imputation(long_data_obj, method = "copyMean", lowerBound = "globalMin", upperBound = "globalMax")
  cat('Imputatron performed copyMeans on your data',"\n")
  Sys.sleep(0.8)
  cat('.')
  Sys.sleep(0.4)
  cat('.')
  Sys.sleep(0.4)
  cat('.',"\n")
  Sys.sleep(0.4)
  cat('###### Imputatron completed data imputation ####',"\n")
  Sys.sleep(0.8)
  cat('Imputatron will start data extraction sequence')
  Sys.sleep(0.4)
  cat('.')
  Sys.sleep(0.4)
  cat('.')
  Sys.sleep(0.4)
  cat('.',"\n")

  # Extract and return data
  df_list <- list()
  for (i in 1:dim(imputed_data@traj)[3]) {
    matrix_data <- imputed_data@traj[, , i]
    df_list[[i]] <- as.data.frame(matrix_data)
  }
  names(df_list) <- imputed_data@varNames
  Sys.sleep(0.8)
  cat('Imputatron extracted a list of dataframes from longData3d object',"\n")

  df_long <- lapply(df_list, function(x) {
    x %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = id_col) %>%
      tidyr::pivot_longer(-dplyr::all_of(id_col), names_to = time_col, values_to = "Value") %>%
      dplyr::mutate(!!rlang::sym(time_col) := as.numeric(gsub("t", "", !!rlang::sym(time_col)))) %>%
      dplyr::arrange(!!rlang::sym(id_col), !!rlang::sym(time_col))
  })
  Sys.sleep(0.8)
  cat('Imputatron melted each dataframe back to long format',"\n")

  df_long <- lapply(seq_along(df_long), function(i) {
    stats::setNames(df_long[[i]], c(id_col, time_col, names(df_list)[i]))
  })
  Sys.sleep(0.8)
  cat('Imputatron added names to the list of dataframes',"\n")
  Sys.sleep(0.8)
  cat('.')
  Sys.sleep(0.4)
  cat('.')
  Sys.sleep(0.4)
  cat('.',"\n")
  Sys.sleep(0.4)
  cat('###### Imputatron completed data extraction ####',"\n")
  Sys.sleep(0.8)
  cat('Imputatron will start tidy dataframe creation sequence')
  Sys.sleep(0.4)
  cat('.')
  Sys.sleep(0.4)
  cat('.')
  Sys.sleep(0.4)
  cat('.',"\n")

  # Combine all data frames
  imputed_df <- Reduce(function(x, y) merge(x, y, by = c(id_col, time_col), all = TRUE), df_long) %>%
    dplyr::arrange(!!rlang::sym(id_col), !!rlang::sym(time_col))
  Sys.sleep(0.8)
  cat('Imputatron created tidy dataframe with all of your variables',"\n")
  Sys.sleep(0.8)
  cat('.')
  Sys.sleep(0.4)
  cat('.')
  Sys.sleep(0.4)
  cat('.',"\n")
  Sys.sleep(0.8)
  cat('Imputatron finished the imputation process.')
  Sys.sleep(0.4)
  cat(paste0('',"\n",
             'You can access the following elements:',"\n",
             "   Imputed Dataframe: $imputed_df", "\n",
             "   Long Data Object: $long_data_obj ", "\n",
             "   Long Data Object with imputed data: $imputed_data ", "\n",
             "   Your original data, but filtered: $filtered_data ", "\n",
             "   Filtered data in wide format: $wide_data ", "\n",
             "   The list with variable names and indexes to create the Long Data Object: $timeInData_list", "\n"))
  Sys.sleep(0.8)
  cat('.')
  Sys.sleep(0.4)
  cat('.')
  Sys.sleep(0.4)
  cat('.',"\n")
  Sys.sleep(0.8)
  cat('Imputatron turning off...',"\n")
  Sys.sleep(0.8)
  cat('~~ Beep Boop Beep Boop ~~ Bye!')


  return(list(imputed_data = imputed_data, imputed_df = imputed_df, filtered_data = filtered_data, wide_data = wide_data, timeInData_list = timeInData_list, long_data_obj = long_data_obj))
}
