#' Run dynamic modelling I-ARIMAX algorithm.
#'
#' @export
#' @importFrom utils globalVariables
#'
#' @param dataframe Your dataframe.
#' @param a_series The first part of the loop.
#' @param b_series The second part of the loop.
#' @param c_series The third part of the loop.
#' @param min_n_subject The minimum number of non NA cases to run the analyses. It will filter cases with more NA's than the threshold. Defaults to 20.
#' @param minvar The minimum variance for both series (&) to include a case. Defaults to 0.01.
#' @param id_var A string containing your id variable.
#' @param metaanalysis Bool to run a random effects meta-analysis or not.
#' @param hlm_compare Optional, to create a comparison with an HLM model, default is FALSE.
#' @param timevar If hlm_compare is TRUE, then a time variable is needed, default is NULL.
#' @param weight_rma If adding an exogenous weight variable to the RMA model.
#' @param weight_rma_var Select weight RMA variable. Defaults as NULL, if NULL then the number of valid observations for Y AND X will be used (!is.na).
#' @param correlation_method Select method for raw semi-partial correlations. Options are: 'spearman', 'pearson' or 'kendall'. Defaults to 'pearson'.
#'
#' @returns a list with dataframes from original iarimax algorithm outuputs + looping machine dataframe.


##############################################
############ Looping Machine function #######
############################################

looping_machine <- function(dataframe,a_series,b_series,c_series, id_var,
                            min_n_subject = 20, minvar = 0.01, hlm_compare = TRUE,
                            metaanalysis = TRUE, timevar = NULL, weight_rma = FALSE,
                            weight_rma_var = NULL, correlation_method = 'pearson') {


  #Create name tags.
  ab_name <- paste0(a_series,"_",b_series)
  bc_name <- paste0(b_series,"_",c_series)
  ca_name <- paste0(c_series,"_",a_series)

  Sys.sleep(1)
  cat(paste0('Calculating a to b: From ',a_series, ' to ',b_series, '...',"\n"))
  Sys.sleep(2)


  #a to b condition.
  a_to_b  <- IARIMAXoid_Pro(dataframe,min_n_subject = min_n_subject, minvar = minvar, y_series = b_series,
                            x_series = a_series, id_var = id_var, hlm_compare = hlm_compare, timevar = timevar)

  a_to_b_sub <- a_to_b$results_df[,c(1,28:31)] #extract.
  a_to_b_sub <- a_to_b_sub %>% #P value.
    dplyr::mutate(p_value = 2 * stats::pt(-abs(xreg / stderr_xreg), n_valid - n_params))
  colnames(a_to_b_sub) <- c(id_var,ab_name,paste0("stderr_",ab_name),paste0(ab_name,"_n_valid"),paste0(ab_name,"_n_params"),paste0(ab_name,'_pval')) #Rename.

  Sys.sleep(1)
  cat(paste0('Calculating b to c: From ' ,b_series, ' to ',c_series, '...',"\n"))
  Sys.sleep(2)

  #b to c condition.
  b_to_c  <- IARIMAXoid_Pro(dataframe,min_n_subject = min_n_subject, minvar = minvar, y_series = c_series,
                            x_series = b_series, id_var = id_var, hlm_compare = hlm_compare, timevar = timevar)

  b_to_c_sub <- b_to_c$results_df[,c(1,28:31)] #extract.
  b_to_c_sub <-   b_to_c_sub %>% #P value.
    dplyr::mutate(p_value = 2 * stats::pt(-abs(xreg / stderr_xreg), n_valid - n_params))
  colnames(b_to_c_sub) <- c(id_var,bc_name,paste0("stderr_",bc_name),paste0(bc_name,"_n_valid"),paste0(bc_name,"_n_params"),paste0(bc_name,'_pval')) #Rename.

  Sys.sleep(1)
  cat(paste0('Calculating c to a: From ',c_series, ' to ',a_series, '...',"\n"))
  Sys.sleep(2)

  #c to a condition.
  c_to_a  <- IARIMAXoid_Pro(dataframe,min_n_subject = min_n_subject, minvar = minvar, y_series = a_series,
                            x_series = c_series, id_var = id_var, hlm_compare = hlm_compare, timevar = timevar)
  c_to_a_sub <- c_to_a$results_df[,c(1,28:31)] #extract.
  c_to_a_sub <-   c_to_a_sub %>% #P value.
    dplyr::mutate(p_value = 2 * stats::pt(-abs(xreg / stderr_xreg), n_valid - n_params))
  colnames(c_to_a_sub) <- c(id_var,ca_name,paste0("stderr_",ca_name),paste0(ca_name,"_n_valid"),paste0(ca_name,"_n_params"),paste0(ca_name,'_pval')) #Rename.


  loop_df <- merge(a_to_b_sub,b_to_c_sub, by = id_var)
  loop_df <- merge(loop_df, c_to_a_sub, by = id_var)

  #Create the loop.
  loop_df$Loop_05_Positive_directed <- ifelse(loop_df[[paste0(ab_name,'_pval')]] < 0.05 & loop_df[[paste0(bc_name,'_pval')]] <0.05 & loop_df[[paste0(ca_name,'_pval')]] <0.05 &
                         loop_df[[ab_name]] > 0 & loop_df[[bc_name]] > 0 & loop_df[[ca_name]] > 0,1,0)

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
  cat(paste('looping Machine Finished.',"\n"))
  cat(paste('',"\n"))
  cat(paste('.'))
  Sys.sleep(0.4)
  cat(paste('.'))
  Sys.sleep(0.4)
  cat(paste('.'))
  Sys.sleep(0.4)
  cat(paste('',"\n"))
  cat(paste('Number of cases with the loop present were:',sum(looping$loop_df$Loop_05_Positive_directed)))


  return(list(loop_df = loop_df,
              iarimax_a_to_b = a_to_b,
              iarimax_b_to_c = b_to_c,
              iarimax_c_to_a = c_to_a))


}

utils::globalVariables(c("xreg", "stderr_xreg", "n_valid", "n_params","ab_name","bc_name","ca_name"))#Declare symbolic global variables.
