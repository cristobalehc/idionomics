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
#' @param hlm_compare Optional, to create a comparison with an HLM model, default is FALSE.
#' @param timevar If hlm_compare is TRUE, then a time variable is needed, default is NULL.
#' @param weight_rma If adding an exogenous weight variable to the RMA model.
#' @param weight_rma_var Select weight RMA variable. Defaults as NULL, if NULL then the number of valid observations for Y AND X will be used (!is.na).
#' @param correlation_method Select method fr raw correlations. Options are: 'spearman', 'pearson' or 'kendall'. Defaults to 'pearson'.
#'
#' @returns A list containing a dataframe with the ARIMA parameters, plus the xreg parameter (the beta value for your x_series) together with their std.errors. If metaanalysis = TRUE, will also output a random effects meta analysis. If hlm_compare = TRUE, will also output a model comparison with HLM.


#######################################
############ I ARIMAX FUNCTION #######
#####################################

IARIMAXoid_Pro <- function(dataframe, min_n_subject = 20, minvar = 0.01, y_series, x_series, id_var,
                           metaanalysis = TRUE, hlm_compare = FALSE, timevar = NULL, weight_rma = FALSE, weight_rma_var = NULL, correlation_method = 'pearson') {



  # CHeck wether variables are in the in the dataset.
  required_vars <- c(y_series, x_series, id_var)

  if (!all(required_vars %in% colnames(dataframe))) {
    missing_vars <- required_vars[!required_vars %in% colnames(dataframe)]
    stop(paste("Cannot find required variables. Check if you spelled the following variables correctly:", paste(missing_vars, collapse = ", ")))
  }

  if (hlm_compare == TRUE) {
    if (is.null(timevar)) {
      stop('You selected hlm_compare, however I cannot compute the model without a time variable. Add it with timevar = "yourtimevariable"')
    }

    if (!(timevar %in% colnames(dataframe))) {
      stop(paste0("The time variable '", timevar, "' was not found in the dataframe. Did you spell it correctly?"))
    }
  }

  if (weight_rma == TRUE & !is.null(weight_rma_var)) {

    if (!(weight_rma_var %in% colnames(dataframe))){

    stop(paste0("The weight variable for RMA '", weight_rma_var, "' was not found in the dataframe. Did you spell it correctly?"))
    }
  }


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
  AR5 <- list()
  stderr_AR5 <- list()

  #MAs
  MA1 <- list()
  stderr_MA1 <- list()
  MA2 <- list()
  stderr_MA2 <- list()
  MA3 <- list()
  stderr_MA3 <- list()
  MA4 <- list()
  stderr_MA4 <- list()
  MA5 <- list()
  stderr_MA5 <- list()

  #Xreg.
  xreg <- list()
  stderr_xreg <- list()

  #Drift.
  drift <- list()
  stderr_drift <- list()

  #Intercept.
  intercept <- list()
  stderr_intercept <- list()

  #Number of valid cases & parameters.
  n_valid <- list()
  n_params <- list()

  #Exclude cases where arimax don't work.
  exclude <- list()

  #Correlations.
  raw_correlation <- list()

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

    #Count number of valid observations.
    n_valid_val <- sum(!is.na(y_vector) & !is.na(x_vector))

    options(warn = -1) #Supress spearman's warning about p-values with ties.
    correlation <- tryCatch(
      {
        stats::cor.test(x = x_vector, y = y_vector, method = correlation_method)
      },
      error = function(e) {
        cat("\n","\n",' Error computing correlation for case: ',as.character(i), "\n","  ",e$message,"\n")
        NULL #Set model as NULL
      }
    )
    options(warn = 0) #Go back to origin.

    #Handle when correlations are null.
    if (is.null(correlation)) {
      raw_correlation[[i]] <- NA
    }


    #Run model and catch errors: TryCatch will do that.
    model <- tryCatch(
      {
        forecast::auto.arima(y = y_vector, xreg = x_vector, approximation = FALSE, stepwise = FALSE)
      },
      error = function(e) {
        cat("\n","\n",'  Error running ARIMAX model for case: ',as.character(i), "\n","  ",e$message,"\n")
        NULL #Set model as NULL
      }
    )

    #Fill the list with NA if the model is null

    if (is.null(model)) {
      AR_N[[i]] <- NA
      I_N[[i]] <- NA
      MA_N[[i]] <- NA
      AR1[[i]] <- NA
      stderr_AR1[[i]] <- NA
      AR2[[i]] <- NA
      stderr_AR2[[i]] <- NA
      AR3[[i]] <- NA
      stderr_AR3[[i]] <- NA
      AR4[[i]] <- NA
      stderr_AR4[[i]] <- NA
      AR5[[i]] <- NA
      stderr_AR5[[i]] <- NA
      MA1[[i]] <- NA
      stderr_MA1[[i]] <- NA
      MA2[[i]] <- NA
      stderr_MA2[[i]] <- NA
      MA3[[i]] <- NA
      stderr_MA3[[i]] <- NA
      MA4[[i]] <- NA
      stderr_MA4[[i]] <- NA
      MA5[[i]] <- NA
      stderr_MA5[[i]] <- NA
      xreg[[i]] <- NA
      stderr_xreg[[i]] <- NA
      drift[[i]] <- NA
      stderr_drift[[i]] <- NA
      intercept[[i]] <- NA
      stderr_intercept[[i]] <- NA

      n_valid[[i]] <- NA
      n_params[[i]] <- NA
      exclude[[i]] <- i
      cat("\n","   Skipping case due to error: ")
      cat("        ... ",round((casen/(length(names))*100),digits = 1),'% completed',"\n","\n") #Keep printing advance percentage.
      next #Skip the next part of this iteration of the loop, so it doesn't get overriden and throws an error.
    }

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
    }    else {
      AR1[[i]] <- NA
      stderr_AR1[[i]] <- NA
    }



    #Fill AR2 parameters conditionally to their existence.
    if ('estimate_ar2' %in% colnames(tidymodel)) {
      AR2[[i]] <- tidymodel$estimate_ar2
      stderr_AR2[[i]] <- tidymodel$std.error_ar2
    }    else {
      AR2[[i]] <- NA
      stderr_AR2[[i]] <- NA
    }


    #Fill AR3 parameters conditionally to their existence.
    if ('estimate_ar3' %in% colnames(tidymodel)) {
      AR3[[i]] <- tidymodel$estimate_ar3
      stderr_AR3[[i]] <- tidymodel$std.error_ar3
    }    else {
      AR3[[i]] <- NA
      stderr_AR3[[i]] <- NA
    }


    #Fill AR4 parameters conditionally to their existence.
    if ('estimate_ar4' %in% colnames(tidymodel)) {
      AR4[[i]] <- tidymodel$estimate_ar4
      stderr_AR4[[i]] <- tidymodel$std.error_ar4
    }     else {
      AR4[[i]] <- NA
      stderr_AR4[[i]] <- NA
    }

    #Fill AR5 parameters conditionally to their existence.
    if ('estimate_ar5' %in% colnames(tidymodel)) {
      AR5[[i]] <- tidymodel$estimate_ar5
      stderr_AR5[[i]] <- tidymodel$std.error_ar5
    }     else {
      AR5[[i]] <- NA
      stderr_AR5[[i]] <- NA
    }


    ##############################
    ###### FIL MA PARAMETERS #####
    #############################

    #Fill MA1 parameters conditionally to their existence.
    if ('estimate_ma1' %in% colnames(tidymodel)) {
      MA1[[i]] <- tidymodel$estimate_ma1
      stderr_MA1[[i]] <- tidymodel$std.error_ma1
    }     else {
      MA1[[i]] <- NA
      stderr_MA1[[i]] <- NA
    }


    #Fill MA2 parameters conditionally to their existence.
    if ('estimate_ma2' %in% colnames(tidymodel)) {
      MA2[[i]] <- tidymodel$estimate_ma2
      stderr_MA2[[i]] <- tidymodel$std.error_ma2
    }     else {
      MA2[[i]] <- NA
      stderr_MA2[[i]] <- NA
    }


    #Fill MA3 parameters conditionally to their existence.
    if ('estimate_ma3' %in% colnames(tidymodel)) {
      MA3[[i]] <- tidymodel$estimate_ma3
      stderr_MA3[[i]] <- tidymodel$std.error_ma3
    }     else {
      MA3[[i]] <- NA
      stderr_MA3[[i]] <- NA
    }


    #Fill MA4 parameters conditionally to their existence.
    if ('estimate_ma4' %in% colnames(tidymodel)) {
      MA4[[i]] <- tidymodel$estimate_ma4
      stderr_MA4[[i]] <- tidymodel$std.error_ma4
    }     else {
      MA4[[i]] <- NA
      stderr_MA4[[i]] <- NA
    }

    #Fill MA4 parameters conditionally to their existence.
    if ('estimate_ma5' %in% colnames(tidymodel)) {
      MA5[[i]] <- tidymodel$estimate_ma5
      stderr_MA5[[i]] <- tidymodel$std.error_ma5
    }     else {
      MA5[[i]] <- NA
      stderr_MA5[[i]] <- NA
    }

    #####################################################
    ##### FILL Xreg, Intercept &  Drift Parameters #####
    ###################################################

    #Fill XREG parameters conditionally to their existence.
    if ('estimate_xreg' %in% colnames(tidymodel)) {
      xreg[[i]] <- tidymodel$estimate_xreg
      stderr_xreg[[i]] <- tidymodel$std.error_xreg
    }     else {
      xreg[[i]] <- NA
      stderr_xreg[[i]] <- NA
    }

    #Fill drift parameters conditionally to their existence.
    if ('estimate_drift' %in% colnames(tidymodel)) {
      drift[[i]] <- tidymodel$estimate_drift
      stderr_drift[[i]] <- tidymodel$std.error_drift
    }     else {
      drift[[i]] <- NA
      stderr_drift[[i]] <- NA
    }

    #Fill intercept parameters conditionally to their existence.
    if ('estimate_intercept' %in% colnames(tidymodel)) {
      intercept[[i]] <- tidymodel$estimate_intercept
      stderr_intercept[[i]] <- tidymodel$std.error_intercept
    }     else {
      intercept[[i]] <- NA
      stderr_intercept[[i]] <- NA
    }

    #Add number of valid cases.
    n_valid[[i]] <- n_valid_val
    n_params[[i]] <- length(model$coef)

    #Fill raw correlation if not null.
    if (!is.null(correlation)) {
      raw_correlation[[i]] <- correlation$estimate[[1]]
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
  AR5_vector <- unlist(AR5)
  stderr_AR5_vector <- unlist(stderr_AR5)
  MA1_vector <- unlist(MA1)
  stderr_MA1_vector <- unlist(stderr_MA1)
  MA2_vector <- unlist(MA2)
  stderr_MA2_vector <- unlist(stderr_MA2)
  MA3_vector <- unlist(MA3)
  stderr_MA3_vector <- unlist(stderr_MA3)
  MA4_vector <- unlist(MA4)
  stderr_MA4_vector <- unlist(stderr_MA4)
  MA5_vector <- unlist(MA5)
  stderr_MA5_vector <- unlist(stderr_MA5)
  xreg_vector <- unlist(xreg)
  stderr_xreg_vector <- unlist(stderr_xreg)
  drift_vector <- unlist(drift)
  stderr_drift_vector <- unlist(stderr_drift)
  intercept_vector <- unlist(intercept)
  stderr_intercept_vector <- unlist(stderr_intercept)
  n_valid_vector <- unlist(n_valid)
  n_params_vector <- unlist(n_params)
  raw_correlation_vector <- unlist(raw_correlation)

  # Combine into a data frame
  results_df <- data.frame(
    Name = names,
    nAR = AR_vector,
    nI = I_vector,
    nMA = MA_vector,
    intercept = intercept_vector,
    AR1 = AR1_vector,
    stderr_AR1 = stderr_AR1_vector,
    AR2 = AR2_vector,
    stderr_AR2 = stderr_AR2_vector,
    AR3 = AR3_vector,
    stderr_AR3 = stderr_AR3_vector,
    AR4 = AR4_vector,
    stderr_AR4 = stderr_AR4_vector,
    AR5 = AR5_vector,
    stderr_AR5 = stderr_AR5_vector,
    MA1 = MA1_vector,
    stderr_MA1 = stderr_MA1_vector,
    MA2 = MA2_vector,
    stderr_MA2 = stderr_MA2_vector,
    MA3 = MA3_vector,
    stderr_MA3 = stderr_MA3_vector,
    MA4 = MA4_vector,
    stderr_MA4 = stderr_MA4_vector,
    MA5 = MA5_vector,
    stderr_MA5 = stderr_MA5_vector,
    drift = drift_vector,
    stderr_drift  = stderr_drift_vector,
    xreg = xreg_vector,
    stderr_xreg = stderr_xreg_vector,
    n_valid = n_valid_vector,
    n_params = n_params_vector,
    raw_correlation = raw_correlation_vector)

  #Set id variable, as id_var for consistency.
  colnames(results_df)[1] <- id_var

  #Run with hlm compare.
  if (hlm_compare == TRUE){


    Sys.sleep(0.8)
    cat(paste('',"\n"))


    Sys.sleep(1.2)
    cat(paste0('',"\n",
               '            1. SUMMARY OF ARIMA PARAMETERS',"\n",
               ' ',"\n",
               " The proportion of AR of order 1 or more is ", round(prop.table(table(results_df$nAR >= 1))[2], 2), "\n",
               " The proportion of I of order 1 or more is ", round(prop.table(table(results_df$nI >= 1))[2], 2), "\n",
               " The proportion of MA of order 1 or more is ", round(prop.table(table(results_df$nMA >= 1))[2], 2), "\n"))

    Sys.sleep(0.8)
    cat(paste('',"\n"))
    cat(paste('',"\n"))
    cat(paste('Running random-effects meta analysis & HLM model for comparison',"\n"))
    cat(paste('',"\n"))

    # IF WEIGHT OPTION IS TRUE.
    if (weight_rma == TRUE & !is.null(weight_rma_var)) {

      weight_rma_var_sym <- rlang::sym(weight_rma_var)

      #Extract first observation of weight variable.
      dataframe_rma_weight <- dataframe %>%
        dplyr::select(!!id_var_sym,!!weight_rma_var_sym) %>% #Subset only weight variable and id variable.
        tidyr::drop_na() %>% #drop_na to avoid extracting a NA weight.
        dplyr::group_by(!!id_var_sym) %>% #group by id.
        dplyr::slice(1L) #Get the first observation.

      #Explicitly set name to merge.
      colnames(dataframe_rma_weight)[1] <- id_var
      colnames(dataframe_rma_weight)[2] <- "weight_variable"

      #merge results_df
      results_df <- merge(results_df, dataframe_rma_weight, by = id_var, all.x = TRUE)

      #Run random effects meta analysis.
      #Try to conduct the random effect meta analysis.
      meta_analysis <-
        tryCatch(
          {
            metafor::rma(yi = results_df$xreg, sei = results_df$stderr_xreg, weights = results_df$weight_variable)
          },
          error = function(e) {
            cat('Error running RME with weight variable:', e$message, '\n','trying with default weight (n of valid observations): \n')

            NULL #Set as null.
          }
        )
      #Print success.
      if (!is.null(meta_analysis)){
        cat(paste0(' Random Effects Meta-Analysis Ran correctly with your specified weight variable: ','"',weight_rma_var,'"\n'))
        ran_rma <- 1 #marker.
      }
      #Delete weight variable.
      results_df <- results_df %>% dplyr::select(-weight_variable)
      #If RMA doesnt work with the especified weight variable.
      if (is.null(meta_analysis)){

        #Get number of valid cases.
        dataframe_count_valid <- dataframe %>%
          dplyr::select(!!id_var_sym,!!y_series_sym,!!x_series_sym) %>% #Subset only id variable and x and y series.
          tidyr::drop_na() %>% #drop_na to avoid extracting a NA weight.
          dplyr::group_by(!!id_var_sym) %>% #group by id.
          dplyr::summarise(weight_variable_n = dplyr::n()) #get the count of valid cases.

        #Merge different weight variable.
        results_df <- merge(results_df, dataframe_count_valid, by = id_var, all.x = TRUE)

        #Try to run the new meta analysis.
        meta_analysis <-
          tryCatch(
            {
              metafor::rma(yi = results_df$xreg, sei = results_df$stderr_xreg, weights = results_df$weight_variable_n)
            },
            error = function(e) {
              cat('Error running RME with default weight (n of valid observations):', e$message, '\n',
                  'Setting weight_rma = FALSE: \n')
              NULL
            }
          )
        #Print success.
        if (!is.null(meta_analysis)){
          cat(' Random Effects Meta-Analysis Ran correctly with the default weight variable (n valid cases)')
        }

        results_df <- results_df %>% dplyr::select(-weight_variable_n)


      }
      #If meta analysis is still null: Just conduct it normally.
      if (is.null(meta_analysis)){
        #Run random effects meta analysis.
        #Try to conduct the random effect meta analysis.
        meta_analysis <-
          tryCatch(
            {
              metafor::rma(yi = results_df$xreg, sei = results_df$stderr_xreg)
            },
            error = function(e) {
              cat('Error running RME:', e$message, '\n')
              NULL
            }
          )

        #Print success.
        if (!is.null(meta_analysis)){
          cat(' Random Effects Meta-Analysis Ran correctly without exogenous weights')
        }

      }

      #If weight RMA equals TRUE, but weight_rma_var is NULL.
    }  else if (weight_rma == TRUE & is.null(weight_rma_var)){



      #Get number of valid cases.
      dataframe_count_valid <- dataframe %>%
        dplyr::select(!!id_var_sym,!!y_series_sym,!!x_series_sym) %>% #Subset only id variable and x and y series.
        tidyr::drop_na() %>% #drop_na to avoid extracting a NA weight.
        dplyr::group_by(!!id_var_sym) %>% #group by id.
        dplyr::summarise(weight_variable_n = dplyr::n()) #get the count of valid cases.

      #Merge different weight variable.
      results_df <- merge(results_df, dataframe_count_valid, by = id_var, all.x = TRUE)

      #Try to run the new meta analysis.
      meta_analysis <-
        tryCatch(
          {
            metafor::rma(yi = results_df$xreg, sei = results_df$stderr_xreg, weights = results_df$weight_variable_n)
          },
          error = function(e) {
            cat('Error running RME with default weight (n of valid observations):', e$message, '\n',
                'Setting weight_rma = FALSE: \n')
            NULL
          }
        )
      if (!is.null(meta_analysis)){
        cat(' Random Effects Meta-Analysis Ran correctly with the default weight variable (n valid cases)')
      }


      #Delete weight variable.
      results_df <- results_df %>% dplyr::select(-weight_variable_n)

      #If meta analysis is still null: Just conduct it normally.
      if (is.null(meta_analysis)){
        #Run random effects meta analysis.
        #Try to conduct the random effect meta analysis.
        meta_analysis <-
          tryCatch(
            {
              metafor::rma(yi = results_df$xreg, sei = results_df$stderr_xreg)
            },
            error = function(e) {
              cat('Error running RME:', e$message, '\n')
              NULL
            }
          )
        if (!is.null(meta_analysis)){
          cat(' Random Effects Meta-Analysis Ran correctly without exogenous weights')
        }
      }

      #If weights_rma = FALSE.
    }   else {



    #Run random effects meta analysis.
      #Try to conduct the random effect meta analysis.
    meta_analysis <-
      tryCatch(
        {
      metafor::rma(yi = results_df$xreg, sei = results_df$stderr_xreg)
        },
      error = function(e) {
        cat('Error running RME:', e$message, '\n')
        NULL
      }
      )
    if (!is.null(meta_analysis)){
      cat(' Random Effects Meta-Analysis Ran correctly without exogenous weights')
    }
    }

    ### LAST CHECK FOR RANDOM META ANALYSIS ###

    #If meta analysis is null: Stop and return early.
    if (is.null(meta_analysis)) {
      cat('Skipping RME and HLM model due to error in RME. Returning metaanalysis = FALSE model. \n')
      return(list(results_df = results_df, error_arimax_skipped = exclude, type = 'contemporaneous_predictor'))
      #If meta analysis worked (not null): Continue with the rest of the calculations (HLM in this case.)
    }
    else {


    Sys.sleep(2)
    cat(paste('',"\n"))
    cat(paste('',"\n"))
    cat(paste('',"\n"))
    cat(paste0('            2. SUMMARY OF RANDOM EFFECTS META ANALYSIS',"\n"))
    print(summary(meta_analysis))

    ##########################
    ##### RUN HLM MODEL #####
    ########################

    #Filter dataframe to be consistent with I-ARIMAX data.
    dataframehlm <- dataframe %>%
      dplyr::filter(!!id_var_sym %in% names) %>% #Filter cases based on minvar and non NA.
      dplyr::filter(!(!!id_var_sym) %in% exclude) #Filter cases excluded due to arimax not working.


    #Run HLM Model.
    # Construct the fixed effects formula
    fixed_formula <- stats::as.formula(paste(y_series, "~", timevar, "+", x_series))

    # Construct the random effects formula
    random_formula <- stats::as.formula(paste("~ 1 +", x_series, "|", id_var))

    # Run HLM Model with tryCatch.
    hlm_model <- tryCatch(
      {
        nlme::lme(fixed = fixed_formula, random = random_formula,
                  data = dataframehlm, na.action = stats::na.omit, correlation = nlme::corAR1(),
                  control = nlme::lmeControl(opt="optim", msMaxIter = 500, niterEM=100, msMaxEval = 500, returnObject = TRUE, tolerance = 1e-06))
      },
      error = function(e) {
        cat('Error running HLM model:', e$message, "\n")
        NULL # Set hlm_model as NULL in case of error
      }
    )

    #Stop early if hlm model had an error.
    if (is.null (hlm_model)){
      # Simplify fixed effects, take away time var.
      fixed_formula <- stats::as.formula(paste(y_series, "~", x_series))
      cat('Simplifying HLM model, deleting time variable from fixed effects',"\n")
      #Run the model again.
      hlm_model <- tryCatch(
        {
          nlme::lme(fixed = fixed_formula, random = random_formula,
                    data = dataframehlm, na.action = stats::na.omit, correlation = nlme::corAR1(),
                    control = nlme::lmeControl(opt="optim", msMaxIter = 500, niterEM=100, msMaxEval = 500, returnObject = TRUE, tolerance = 1e-06))
        },
        error = function(e) {
          cat('Error running simplified hlm model:', e$message, "\n")
          NULL # Set hlm_model as NULL in case of error
        }
      )

    }

    #Stop early if hlm model had an error.
    if (is.null(hlm_model)){
      cat('Skipping HLM model calculations due to error. Returning hlm_model = FALSE model. \n')
      return(list(results_df = results_df,meta_analysis = meta_analysis, error_arimax_skipped = exclude, type = 'contemporaneous_predictor'))
      #If hlm model is not null, then do the rest of the calculations.
    }
    else {



      # Create dataframe with random effects.
      df_rand <- as.data.frame(hlm_model$coefficients$random) #Create dataframe
      df_rand <- tibble::rownames_to_column(df_rand) #Add rowname as column
      colnames(df_rand) <- c(id_var,'Intercept',x_series) #Change names to legible ones.
      df_rand$random_slope <- df_rand[[x_series]]+hlm_model$coefficients$fixed[3] #Create random slopes.

      #Values for comparison.
      iarimaxest <- meta_analysis$b[1] #Extract random meta analysis estimate.
      iarimaxse <- meta_analysis$se #Extract random meta analysis se.
      hlmest <- hlm_model$coefficients$fixed[[3]] #Extract the fixed effect of x_series, the order depends on the formula.
      hlmse <- sqrt(hlm_model$varFix[9]) #Extract the se, for the fixed effect of x_series, the order depends on the formula.
      iarimaxmean <- mean(results_df$xreg, na.rm = TRUE) #Calculate mean of individual xregs.
      iarimaxvar <- stats::var(results_df$xreg, na.rm = TRUE) #Calculate variance of individual xregs.
      corrmean <- mean(results_df$raw_correlation, na.rm = TRUE) #Calculate mean of raw correlations.
      corrvar <- stats::var(results_df$raw_correlation, na.rm = TRUE) #Calculate variance of raw correlations.
      hlmmean <- mean(df_rand$random_slope, na.rm = TRUE) #Calculate mean of random slopes.
      hlmvar <- stats::var(df_rand$random_slope, na.rm = TRUE) #Calculate variance of random slopes.

      #Opposite cases.
      num_oppcase_iarimax <- ifelse(meta_analysis$b[1] > 0,sum(results_df$xreg <0, na.rm = TRUE),sum(results_df$xreg > 0, na.rm = TRUE))
      num_oppcase_hlm <- ifelse(hlm_model$coefficients$fixed[[3]] > 0,sum(df_rand$random_slope <0, na.rm = TRUE),sum(df_rand$random_slope > 0, na.rm = TRUE))
      num_oppcase_corr <- ifelse(meta_analysis$b[1] > 0,sum(results_df$raw_correlation <0, na.rm = TRUE),sum(results_df$raw_correlation > 0, na.rm = TRUE))

      iarimaxtohlm <- list(IarimaxEstimate = iarimaxest,
                           IarimaxSE = iarimaxse,
                           HLMEstimate = hlmest,
                           HLMSE = hlmse,
                           IarimaxMean = iarimaxmean,
                           IarimaxVariance = iarimaxvar,
                           HLMMean = hlmmean,
                           HLMVariance = hlmvar,
                           Corrmean = corrmean,
                           Corrvar = corrvar,
                           OppositesIarimax = num_oppcase_iarimax,
                           OppositesHLM = num_oppcase_hlm,
                           OppositesCorr = num_oppcase_corr)
      Sys.sleep(2)
      #Print Summary.
      cat(paste('',"\n"))
      cat(paste('',"\n"))
      cat(paste0('            3. I-ARIMAX V/S HLM COMPARISON ',"\n",
                 '   ',"\n",
                 '   OMNIBUS - RMA ESTIMATES VS HLM',"\n",
                 "\n",
                 '            I-ARIMAX                         HLM     ',"\n",
                 '   ----------------------------------------------------',"\n",
                 '   IARIMAX est.   IARIMAX s.e.       HLM est.   HLM s.e.',"\n",
                 '      ',round(iarimaxest,3),
                 '          ',round(iarimaxse,3),'            ',
                 round(hlmest,3),'      ',round(hlmse,3),"\n",
                 "\n",
                 '   SUMMARY OF INDIVIDUAL EFFECTS',"\n",
                 "\n",
                 '    AVERAGES:',"\n",
                 '                            MEAN:',"\n",
                 '     * IARIMAX mean effect:                       ',round(iarimaxmean,4),"\n",
                 '     * HLM mean effect:                           ',round(hlmmean,4),"\n",
                 '     * Raw corr mean effect:                      ',round(corrmean,4),"\n",
                 '                          VARIANCE:',"\n",
                 '     * IARIMAX variance of individual effects:    ',round(iarimaxvar,4),"\n",
                 '     * HLM variance of random effects:            ',round(hlmvar,4),"\n",
                 '     * Raw corr variance:                         ',round(corrvar,4),"\n","\n",
                 '    OPPOSITES:',"\n",
                 '     * IARIMAX N cases with opposite direction:   ',round(num_oppcase_iarimax,4),"\n",
                 '     * HLM N cases with opposite direction:       ',round(num_oppcase_hlm,4),"\n",
                 '     * Corr N cases with opposite direction:      ',round(num_oppcase_corr,4)))


      Sys.sleep(1.2)
      cat(paste('',"\n"))
      cat(paste('.'))
      Sys.sleep(0.4)
      cat(paste('.'))
      Sys.sleep(0.4)
      cat(paste('.'))
      Sys.sleep(0.4)
      cat(paste('',"\n"))
      Sys.sleep(0.8)
      cat(paste('I-ARIMAX algorithm finished.',"\n"))
      cat(paste('',"\n"))

      #Return values.
      return(list(results_df = results_df,meta_analysis = meta_analysis, hlm_mod = hlm_model, rand_df = df_rand, comparison = iarimaxtohlm, error_arimax_skipped = exclude, type = 'contemporaneous_predictor'))
     }
    }

  }
  #Run just the meta analysis (This should work only if meta analysis was not null, if null it stopped before)
  else if (metaanalysis == TRUE & hlm_compare == FALSE){


    Sys.sleep(0.8)
    cat(paste('',"\n"))
    cat(paste('Running random-effects meta analysis.',"\n"))
    cat(paste('',"\n"))

    # IF WEIGHT OPTION IS TRUE.
    if (weight_rma == TRUE & !is.null(weight_rma_var)) {

      weight_rma_var_sym <- rlang::sym(weight_rma_var)

      #Extract first observation of weight variable.
      dataframe_rma_weight <- dataframe %>%
        dplyr::select(!!id_var_sym,!!weight_rma_var_sym) %>% #Subset only weight variable and id variable.
        tidyr::drop_na() %>% #drop_na to avoid extracting a NA weight.
        dplyr::group_by(!!id_var_sym) %>% #group by id.
        dplyr::slice(1L) #Get the first observation.

      #Explicitly set name to merge.
      colnames(dataframe_rma_weight)[1] <- id_var
      colnames(dataframe_rma_weight)[2] <- "weight_variable"



      #merge results_df
      results_df <- merge(results_df, dataframe_rma_weight, by = id_var, all.x = TRUE)


      #Run random effects meta analysis.
      #Try to conduct the random effect meta analysis.
      meta_analysis <-
        tryCatch(
          {
            metafor::rma(yi = results_df$xreg, sei = results_df$stderr_xreg, weights = results_df$weight_variable)
          },
          error = function(e) {
            cat('Error running RME with weight variable:', e$message, '\n','trying with default weight (n of valid observations): \n')

            NULL #Set as null.
          }
        )
      #Print success.
      if (!is.null(meta_analysis)){
        cat(paste0(' Random Effects Meta-Analysis Ran correctly with your specified weight variable: ','"',weight_rma_var,'"\n'))
        ran_rma <- 1 #marker.
      }
      #Delete weight variable.
      results_df <- results_df %>% dplyr::select(-weight_variable)
      #If RMA doesnt work with the especified weight variable.
      if (is.null(meta_analysis)){

        #Get number of valid cases.
        dataframe_count_valid <- dataframe %>%
          dplyr::select(!!id_var_sym,!!y_series_sym,!!x_series_sym) %>% #Subset only id variable and x and y series.
          tidyr::drop_na() %>% #drop_na to avoid extracting a NA weight.
          dplyr::group_by(!!id_var_sym) %>% #group by id.
          dplyr::summarise(weight_variable_n = dplyr::n()) #get the count of valid cases.

        #Merge different weight variable.
        results_df <- merge(results_df, dataframe_count_valid, by = id_var, all.x = TRUE)

        #Try to run the new meta analysis.
        meta_analysis <-
          tryCatch(
            {
              metafor::rma(yi = results_df$xreg, sei = results_df$stderr_xreg, weights = results_df$weight_variable_n)
            },
            error = function(e) {
              cat('Error running RME with default weight (n of valid observations):', e$message, '\n',
                  'Setting weight_rma = FALSE: \n')
              NULL
            }
          )
        #Print success.
        if (!is.null(meta_analysis)){
          cat(' Random Effects Meta-Analysis Ran correctly with the default weight variable (n valid cases)')
        }

        results_df <- results_df %>% dplyr::select(-weight_variable_n)


      }
      #If meta analysis is still null: Just conduct it normally.
      if (is.null(meta_analysis)){
        #Run random effects meta analysis.
        #Try to conduct the random effect meta analysis.
        meta_analysis <-
          tryCatch(
            {
              metafor::rma(yi = results_df$xreg, sei = results_df$stderr_xreg)
            },
            error = function(e) {
              cat('Error running RME:', e$message, '\n')
              NULL
            }
          )

        #Print success.
        if (!is.null(meta_analysis)){
          cat(' Random Effects Meta-Analysis Ran correctly without exogenous weights')
        }

      }

      #If weight RMA equals TRUE, but weight_rma_var is NULL.
    }  else if (weight_rma == TRUE & is.null(weight_rma_var)){



      #Get number of valid cases.
      dataframe_count_valid <- dataframe %>%
        dplyr::select(!!id_var_sym,!!y_series_sym,!!x_series_sym) %>% #Subset only id variable and x and y series.
        tidyr::drop_na() %>% #drop_na to avoid extracting a NA weight.
        dplyr::group_by(!!id_var_sym) %>% #group by id.
        dplyr::summarise(weight_variable_n = dplyr::n()) #get the count of valid cases.

      #Merge different weight variable.
      results_df <- merge(results_df, dataframe_count_valid, by = id_var, all.x = TRUE)

      #Try to run the new meta analysis.
      meta_analysis <-
        tryCatch(
          {
            metafor::rma(yi = results_df$xreg, sei = results_df$stderr_xreg, weights = results_df$weight_variable_n)
          },
          error = function(e) {
            cat('Error running RME with default weight (n of valid observations):', e$message, '\n',
                'Setting weight_rma = FALSE: \n')
            NULL
          }
        )
      if (!is.null(meta_analysis)){
        cat(' Random Effects Meta-Analysis Ran correctly with the default weight variable (n valid cases)')
      }


      #Delete weight variable.
      results_df <- results_df %>% dplyr::select(-weight_variable_n)

      #If meta analysis is still null: Just conduct it normally.
      if (is.null(meta_analysis)){
        #Run random effects meta analysis.
        #Try to conduct the random effect meta analysis.
        meta_analysis <-
          tryCatch(
            {
              metafor::rma(yi = results_df$xreg, sei = results_df$stderr_xreg)
            },
            error = function(e) {
              cat('Error running RME:', e$message, '\n')
              NULL
            }
          )
        if (!is.null(meta_analysis)){
          cat(' Random Effects Meta-Analysis Ran correctly without exogenous weights')
        }
      }

      #If weights_rma = FALSE.
    }   else {



      #Run random effects meta analysis.
      #Try to conduct the random effect meta analysis.
      meta_analysis <-
        tryCatch(
          {
            metafor::rma(yi = results_df$xreg, sei = results_df$stderr_xreg)
          },
          error = function(e) {
            cat('Error running RME:', e$message, '\n')
            NULL
          }
        )
      if (!is.null(meta_analysis)){
        cat(' Random Effects Meta-Analysis Ran correctly without exogenous weights')
      }
    }
    #####################################
    #LAST CHECK OF RANDOM META ANALYSIS#
    ###################################
    #If meta analysis is null: Stop and return early.
    if (is.null(meta_analysis)) {
      cat('Skipping RME due to error. Returning metaanalysis = FALSE model. \n')
      return(list(results_df = results_df, error_arimax_skipped = exclude, type = 'contemporaneous_predictor'))

    }
    else {

      Sys.sleep(0.8)
      cat(paste('',"\n"))


    cat(paste0('',"\n",
               '            1. SUMMARY OF ARIMA PARAMETERS',"\n",
               ' ',"\n",
               " The proportion of AR of order 1 or more is ", round(prop.table(table(results_df$nAR >= 1))[2], 2), "\n",
               " The proportion of I of order 1 or more is ", round(prop.table(table(results_df$nI >= 1))[2], 2), "\n",
               " The proportion of MA of order 1 or more is ", round(prop.table(table(results_df$nMA >= 1))[2], 2), "\n"))

    Sys.sleep(0.8)
    cat(paste('',"\n"))
    cat(paste0('            2. SUMMARY OF RANDOM EFFECTS META ANALYSIS',"\n"))
    print(summary(meta_analysis))
    Sys.sleep(1.2)
    cat(paste('',"\n"))
    cat(paste('.'))
    Sys.sleep(0.4)
    cat(paste('.'))
    Sys.sleep(0.4)
    cat(paste('.'))
    Sys.sleep(0.4)
    cat(paste('',"\n"))
    Sys.sleep(0.8)
    cat(paste('I-ARIMAX algorithm finished.',"\n"))
    cat(paste('',"\n"))

    return(list(results_df = results_df,meta_analysis = meta_analysis, error_arimax_skipped = exclude, type = 'contemporaneous_predictor'))

    }
  }
  else {


    Sys.sleep(0.8)
    cat(paste('',"\n"))

    cat(paste0('',"\n",
               '            1. SUMMARY OF ARIMA PARAMETERS',"\n",
               "The proportion of AR of order 1 or more is ", round(prop.table(table(results_df$nAR >= 1))[2], 2), "\n",
               "The proportion of I of order 1 or more is ", round(prop.table(table(results_df$nI >= 1))[2], 2), "\n",
               "The proportion of MA of order 1 or more is ", round(prop.table(table(results_df$nMA >= 1))[2], 2), "\n"))

    Sys.sleep(0.8)
    cat(paste('',"\n"))
    cat(paste0('2. No Random-Effects meta-analysis or HLM model comparison were made.',"\n"))
    Sys.sleep(1.2)
    cat(paste('',"\n"))
    cat(paste('.'))
    Sys.sleep(0.4)
    cat(paste('.'))
    Sys.sleep(0.4)
    cat(paste('.'))
    Sys.sleep(0.4)
    cat(paste('',"\n"))
    Sys.sleep(0.8)
    cat(paste('I-ARIMAX algorithm finished.',"\n"))
    cat(paste('',"\n"))
    return(list(results_df = results_df, error_arimax_skipped = exclude, type = 'contemporaneous_predictor'))

  }

}


utils::globalVariables(c("count", "var_y", "var_x",'weight_variable','weight_variable_n')) #Declare symbolic global variables.
