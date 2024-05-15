#' Run I-ARIMAX algorithm.
#'
#' @export
#' @importFrom utils globalVariables
#'
#' @param dataframe Your dataframe.
#' @param min_n_subject The minimum number of non NA cases to run the analyses. It will filter cases with more NA's than the threshold. Defaults to 20.
#' @param minvar The minimum variance for both series (&) to include a case. Defaults to 0.01.
#' @param y_series A string containing the name of your dependent variable y.
#' @param x_series A string containing the name of your independent variable x.
#' @param id_var A string containing your id variable.
#' @param metaanalysis Bool to run a random effects meta-analysis or not.
#'
#' @returns A list containing a dataframe with the ARIMA parameteres, plus the xreg parameter (the beta value for your x_series) together with their std.errors. If metaanalysis = TRUE, will also output a random effects meta analysis.

#######################################
############ I ARIMAX FUNCTION #######
#####################################

IARIMAXoid_Pro <- function(dataframe, min_n_subject = 20, minvar = 0.01, y_series, x_series, id_var, metaanalysis = TRUE) {


    # dataframe = your dataframe's name.
    # min_n_subject = The minimum number of non na cases to run the analyses.
    # minvar = minimum variance to include a case.
    # y_series = the name of your dependent variable, as string.
    # y_series = the name of your independent variable, as string.
    # id_var = the name of your id varible per subject within your dataframe, as string.
    # metaanalysis = bool to run a random effects meta-analysis or not.

    # Convert strings to rlang::symbols
    y_series_sym <- rlang::sym(y_series)
    x_series_sym <- rlang::sym(x_series)
    id_var_sym <- rlang::sym(id_var)

    cat(paste('Filtering your data based on specified minimum non NA per subject and variance'))
    Sys.sleep(0.4)
    cat(paste('.'))
    Sys.sleep(0.4)
    cat(paste('.'))
    Sys.sleep(0.4)
    cat(paste('.',"\n"))
    Sys.sleep(0.8)
    cat(paste('',"\n"))

    # Filter N complete observations with variance conditions
    counts <- dataframe %>%
      dplyr::group_by(!!id_var_sym) %>%
      dplyr::filter(!is.na(!!y_series_sym) & !is.na(!!x_series_sym)) %>%
      dplyr::summarise(
        count = dplyr::n(),
        var_y = stats::var(!!y_series_sym, na.rm = TRUE),
        var_x = stats::var(!!x_series_sym, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      dplyr::filter(count >= min_n_subject, var_y >= minvar, var_x >= minvar)


    # Get the names of the id_var, and force them as characters
    names <- as.character(counts[[id_var]])

    cat(paste('Your data was filtered:',length(names),'subjects will be used for the analyses',"\n"))
    Sys.sleep(0.8)


    # Create lists to terms
    AR_N <- list()
    I_N <- list()
    MA_N <- list()

    #Create lists to store the parameters.

    #ARs.
    AR1 <- list()
    stderr_AR1 <- list()
    AR2 <- list()
    stderr_AR2 <- list()
    AR3 <- list()
    stderr_AR3 <- list()
    AR4 <- list()
    stderr_AR4 <- list()

    #MAs
    MA1 <- list()
    stderr_MA1 <- list()
    MA2 <- list()
    stderr_MA2 <- list()
    MA3 <- list()
    stderr_MA3 <- list()
    MA4 <- list()
    stderr_MA4 <- list()


    #Xreg.
    xreg <- list()
    stderr_xreg <- list()


    # Run auto.arima per case
    cat(paste('',"\n"))
    cat(paste('Running I-ARIMAX algorithm'))
    Sys.sleep(0.4)
    cat(paste('.'))
    Sys.sleep(0.4)
    cat(paste('.'))
    Sys.sleep(0.4)
    cat(paste('.',"\n"))
    Sys.sleep(0.4)
    cat(paste('',"\n"))
    Sys.sleep(0.8)




    #Start case number counter.
    casen = 0
    for(i in names) {

      #Update case number.
      casen <- casen + 1

      #Start text.
      cat(paste('  Applying auto ARIMAX to case: ', as.character(i), ' ... '))

      ##################################################
      ############### DEV COMMENT ######################
      ######## HANDLE NA FOR EAC VECTOR BELOW #########
      #################################################

      #Extract y vector.
      y_vector <- dataframe %>%
        dplyr::filter(!!id_var_sym == i) %>%
        dplyr::pull(!!y_series_sym) #Should I add NA omit? I think so...

      #Extract x vector.
      x_vector <- dataframe %>%
        dplyr::filter(!!id_var_sym == i) %>%
        dplyr::pull(!!x_series_sym) #Should I add NA omit? I think so...

      #Run model.
      model <- forecast::auto.arima(y = y_vector, xreg = x_vector, approximation = FALSE, stepwise = FALSE)
      #Tidy the model.
      tidymodel <- broom::tidy(model)
      #Cast the tidy dataframe.
      tidymodel <-  tidyr::pivot_wider(tidymodel, names_from ='term', values_from=c('estimate','std.error'))

      #Fill the number of ARIMA parameteres lists: Just the number of AR I MA processes involved.
      AR_N[[i]] <- model$arma[1] #Fill AR list.
      I_N[[i]] <- model$arma[6] #Fill I list.
      MA_N[[i]] <- model$arma[2] #Fill MA List.

      ##############################
      ###### FIL AR PARAMETERS #####
      #############################

      #Fill AR1 parameters conditionally to their existence.
      if ('estimate_ar1' %in% colnames(tidymodel)) {
        AR1[[i]] <- tidymodel$estimate_ar1
        stderr_AR1[[i]] <- tidymodel$std.error_ar1
      }
      else {
        AR1[[i]] <- NA
        stderr_AR1[[i]] <- NA
      }


    #Fill AR2 parameters conditionally to their existence.
    if ('estimate_ar2' %in% colnames(tidymodel)) {
      AR2[[i]] <- tidymodel$estimate_ar2
      stderr_AR2[[i]] <- tidymodel$std.error_ar2
    }
    else {
      AR2[[i]] <- NA
      stderr_AR2[[i]] <- NA
    }

    #Fill AR3 parameters conditionally to their existence.
    if ('estimate_ar3' %in% colnames(tidymodel)) {
      AR3[[i]] <- tidymodel$estimate_ar3
      stderr_AR3[[i]] <- tidymodel$std.error_ar3
    }
    else {
      AR3[[i]] <- NA
      stderr_AR3[[i]] <- NA
    }

    #Fill AR4 parameters conditionally to their existence.
    if ('estimate_ar4' %in% colnames(tidymodel)) {
      AR4[[i]] <- tidymodel$estimate_ar4
      stderr_AR4[[i]] <- tidymodel$std.error_ar4
    }
    else {
      AR4[[i]] <- NA
      stderr_AR4[[i]] <- NA
    }

    ##############################
    ###### FIL MA PARAMETERS #####
    #############################

    #Fill MA1 parameters conditionally to their existence.
    if ('estimate_ma1' %in% colnames(tidymodel)) {
      MA1[[i]] <- tidymodel$estimate_ma1
      stderr_MA1[[i]] <- tidymodel$std.error_ma1
    }
    else {
      MA1[[i]] <- NA
      stderr_MA1[[i]] <- NA
    }

    #Fill MA2 parameters conditionally to their existence.
    if ('estimate_ma2' %in% colnames(tidymodel)) {
      MA2[[i]] <- tidymodel$estimate_ma2
      stderr_MA2[[i]] <- tidymodel$std.error_ma2
    }
    else {
      MA2[[i]] <- NA
      stderr_MA2[[i]] <- NA
    }

    #Fill MA3 parameters conditionally to their existence.
    if ('estimate_ma3' %in% colnames(tidymodel)) {
      MA3[[i]] <- tidymodel$estimate_ma3
      stderr_MA3[[i]] <- tidymodel$std.error_ma3
    }
    else {
      MA3[[i]] <- NA
      stderr_MA3[[i]] <- NA
    }

    #Fill MA4 parameters conditionally to their existence.
    if ('estimate_ma4' %in% colnames(tidymodel)) {
      MA4[[i]] <- tidymodel$estimate_ma4
      stderr_MA4[[i]] <- tidymodel$std.error_ma4
    }
    else {
      MA4[[i]] <- NA
      stderr_MA4[[i]] <- NA
    }

    #Fill XREG parameters conditionally to their existence.
    if ('estimate_xreg' %in% colnames(tidymodel)) {
      xreg[[i]] <- tidymodel$estimate_xreg
      stderr_xreg[[i]] <- tidymodel$std.error_xreg
    }
    else {
      xreg[[i]] <- NA
      stderr_xreg[[i]] <- NA
    }


    #Finish the text.
    cat(round((casen/(length(names))*100),digits = 1),'% completed',"\n")
  }


  ###############################################
  ######## CREATE DATAFRAME TO RETURN ###########
  ###############################################

    # Convert the lists to vectors
    AR_vector <- unlist(AR_N)
    I_vector <- unlist(I_N)
    MA_vector <- unlist(MA_N)
    AR1_vector <- unlist(AR1)
    stderr_AR1_vector <- unlist(stderr_AR1)
    AR2_vector <- unlist(AR2)
    stderr_AR2_vector <- unlist(stderr_AR2)
    AR3_vector <- unlist(AR3)
    stderr_AR3_vector <- unlist(stderr_AR3)
    AR4_vector <- unlist(AR4)
    stderr_AR4_vector <- unlist(stderr_AR4)
    MA1_vector <- unlist(MA1)
    stderr_MA1_vector <- unlist(stderr_MA1)
    MA2_vector <- unlist(MA2)
    stderr_MA2_vector <- unlist(stderr_MA2)
    MA3_vector <- unlist(MA3)
    stderr_MA3_vector <- unlist(stderr_MA3)
    MA4_vector <- unlist(MA4)
    stderr_MA4_vector <- unlist(stderr_MA4)
    xreg_vector <- unlist(xreg)
    stderr_xreg_vector <- unlist(stderr_xreg)


    # Combine into a data frame
    results_df <- data.frame(
      Name = names,
      nAR = AR_vector,
      nI = I_vector,
      nMA = MA_vector,
      AR1 = AR1_vector,
      stderr_AR1 = stderr_AR1_vector,
      AR2 = AR2_vector,
      stderr_AR2 = stderr_AR2_vector,
      AR3 = AR3_vector,
      stderr_AR3 = stderr_AR3_vector,
      AR4 = AR4_vector,
      stderr_AR4 = stderr_AR4_vector,
      MA1 = MA1_vector,
      stderr_MA1 = stderr_MA1_vector,
      MA2 = MA2_vector,
      stderr_MA2 = stderr_MA2_vector,
      MA3 = MA3_vector,
      stderr_MA3 = stderr_MA3_vector,
      MA4 = MA4_vector,
      stderr_MA4 = stderr_MA4_vector,
      xreg = xreg_vector,
      stderr_xreg = stderr_xreg_vector)

    if (metaanalysis == TRUE){
    Sys.sleep(0.8)
    cat(paste('',"\n"))
    cat(paste('Running random-effects meta analysis.',"\n"))
    cat(paste('',"\n"))

    meta_analysis <- metafor::rma(yi = results_df$xreg, sei = results_df$stderr_xreg)

    Sys.sleep(0.8)
    cat(paste('I-ARIMAX algorithm finished.',"\n"))
    cat(paste('',"\n"))


    cat(paste0('',"\n",
               '1. Summary of ARIMA parameters.',"\n",
               ' ',"\n",
               " The proportion of AR of order 1 or more is ", round(prop.table(table(results_df$nAR >= 1))[2], 2), "\n",
               " The proportion of I of order 1 or more is ", round(prop.table(table(results_df$nI >= 1))[2], 2), "\n",
               " The proportion of MA of order 1 or more is ", round(prop.table(table(results_df$nMA >= 1))[2], 2), "\n"))

    Sys.sleep(0.8)
    cat(paste('',"\n"))
    cat(paste0('2. Summary of Random-Effects meta-analysis.',"\n"))
    print(summary(meta_analysis))

    return(list(results_df = results_df,meta_analysis = meta_analysis))
    }
    else {

      Sys.sleep(0.8)
      cat(paste('I-ARIMAX algorithm finished.',"\n"))
      cat(paste('',"\n"))

      cat(paste0('',"\n",
                 '1. Summary of ARIMA parameters:',"\n",
                 "The proportion of AR of order 1 or more is ", round(prop.table(table(results_df$nAR >= 1))[2], 2), "\n",
                 "The proportion of I of order 1 or more is ", round(prop.table(table(results_df$nI >= 1))[2], 2), "\n",
                 "The proportion of MA of order 1 or more is ", round(prop.table(table(results_df$nMA >= 1))[2], 2), "\n"))

      Sys.sleep(0.8)
      cat(paste('',"\n"))
      cat(paste0('2. No Random-Effects meta-analysis was done.',"\n"))
    return(list(results_df = results_df))

    }
  }

utils::globalVariables(c("count", "var_y", "var_x")) #Declare symbolic global variables.

