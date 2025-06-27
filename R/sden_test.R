#' Run Sign Divergence or Equisyncratic Null tests.
#'
#' @export
#'
#' @param iarimax_object Your iarimax object.
#' @param alpha_arimax Critical value for arimax parameter significance. Defaults to 0.05.
#' @param alpha_binom Critical value for binomial test. Defaults to alpha_arimax for ENT, and alpha_arimax/2 for SDT
#' @param test type of test, default to "auto". "auto", "SDT", and "ENT" are supported.
#'
#' @returns The appropriate sden test.


#######################################
############ sden test     ############
#####################################



sden_test <- function(iarimax_object,alpha_arimax = 0.05, alpha_binom = NULL, test = "auto"){


  #Guard clause.

  # Check if test is valid
  if (!test %in% c('auto', 'SDT', 'ENT')) {
    stop('Invalid test type. Supported tests are "auto", "SDT", and "ENT".')
  }

  # Check if arimax_object is an iarimax object
  if (!is.list(iarimax_object) || !all(c('results_df', 'meta_analysis') %in% names(iarimax_object))) {
    stop('sden_test requieres an iarimax_object object that contains "results_df" and "meta_analysis". Run IARIMAXoid_Pro with metaanalysis = TRUE')
  }

  #Start the test.

  # Get p-values based on t-distribution.
  iarimax_object$results_df$pval <- 2*stats::pt(-abs(iarimax_object$results_df$xreg /
                                                       iarimax_object$results_df$stderr_xreg),
                                                iarimax_object$results_df$n_valid - iarimax_object$results_df$n_params)

  # Get number of valid effects (with p-values)
  number_of_effects <- sum(!is.na(iarimax_object$results_df$pval))

  # Get REMA parameters.
  rema_beta <- iarimax_object$meta_analysis$beta
  rema_pval <- iarimax_object$meta_analysis$pval

  # Get number of statistically signficant cases.
  positive_sig_sum <- sum(ifelse(iarimax_object$results_df$xreg > 0 &
                                   iarimax_object$results_df$pval < alpha_arimax,1,0), na.rm = TRUE)
  negative_sig_sum <- sum(ifelse(iarimax_object$results_df$xreg < 0 &
                                   iarimax_object$results_df$pval < alpha_arimax,1,0), na.rm = TRUE)
  all_sig_sum <- sum(ifelse(iarimax_object$results_df$pval < alpha_arimax,1,0), na.rm = TRUE)


  #Automatically select the test.
  if (test == "auto") {

    if (rema_pval >= 0.05) {
      message('Running Equisyncratic Null Test (ENT)')

      sig_effects <-  all_sig_sum
      pnull <- ifelse(is.null(alpha_binom), alpha_arimax, alpha_binom)
      ttype <- 'ENT'

    }
    else if (rema_pval < 0.05 & rema_beta > 0) {
      message('Running Sign Divergence Test (SDT): number of negative cases (counter positive pooled-effect) values are being evaluated')

      sig_effects <- negative_sig_sum
      pnull <- ifelse(is.null(alpha_binom), alpha_arimax/2, alpha_binom)
      ttype <- 'SDT counter-positive'

    }
    else if (rema_pval < 0.05 & rema_beta < 0) {

      message('Running Sign Divergence Test (SDT): number of positive cases (counter negative pooled-effect) values are being evaluated')

      sig_effects <- positive_sig_sum
      pnull <- ifelse(is.null(alpha_binom), alpha_arimax/2, alpha_binom)
      ttype <- 'SDT counter-negative'

    }
    else {
      message('Meta-Analysis beta value is exactly zero, defaulting to ENT')
      message('Running Equisyncratic Null Test (ENT)')

      if (rema_pval < 0.05) {
        message(paste0('The p-value of the pooled effect statistically significant at 0.05.\n p = ', round(rema_pval,4), ".\n SDT test is suggested"))
      }

      sig_effects <- all_sig_sum
      pnull <- ifelse(is.null(alpha_binom), alpha_arimax, alpha_binom)
      ttype <- 'ENT'

    }

  }

  else if (test == "SDT") {

    if (rema_beta > 0) {
      message('Running Sign Divergence Test (SDT): number of negative cases (counter positive pooled-effect) values are being evaluated')

      if (rema_pval >= 0.05) {
        message(paste0('The p-value of the pooled effect is not significant at 0.05.\n p = ', round(rema_pval,4), ". \n ENT test is suggested"))
      }

      sig_effects <- negative_sig_sum
      pnull <- ifelse(is.null(alpha_binom), alpha_arimax/2, alpha_binom)
      ttype <- 'SDT counter-positive'


    }
    else if (rema_beta < 0) {

      if (rema_pval >= 0.05) {
        message(paste0('The p-value of the pooled effect is not significant at 0.05.\n p = ', round(rema_pval,4), ". \n ENT test is suggested"))
      }

      message('Running Sign Divergence Test (SDT): number of positive cases (counter negative pooled-effect) values are being evaluated')

      sig_effects <- positive_sig_sum
      pnull <- ifelse(is.null(alpha_binom), alpha_arimax/2, alpha_binom)
      ttype <- 'SDT counter-negative'


    }
    else {
      message('Meta-Analysis beta value is exactly zero, defaulting to ENT')
      message('Running Equisyncratic Null Test (ENT)')

      if (rema_pval < 0.05) {
        message(paste0('The p-value of the pooled effect statistically significant at 0.05.\n p = ', round(rema_pval,4), ".\n SDT test is suggested"))
      }

      sig_effects <- all_sig_sum
      pnull <- ifelse(is.null(alpha_binom), alpha_arimax, alpha_binom)
      ttype <- 'ENT'

    }

  }
  else if (test == "ENT") {

    message('Running Equisyncratic Null Test (ENT)')

    if (rema_pval < 0.05) {
      message(paste0('The p-value of the pooled effect statistically significant at 0.05.\n p = ', round(rema_pval,4), ".\n SDT test is suggested"))
    }

    sig_effects <- all_sig_sum
    pnull <- ifelse(is.null(alpha_binom), alpha_arimax, alpha_binom)
    ttype <- 'ENT'

  }
  #Run the appropriate binomial test.

  bintest <- stats::binom.test(x = sig_effects, n = number_of_effects, p = pnull, alternative = "greater")

  #Structure and present output.

  cat(paste0("\n","\n",
             '                  SDEN Test Results ',"\n",
             '               -----------------------',"\n",
             '   Test: ',ttype, ifelse(ttype == 'ENT',' (Equisyncratic Null Test)',' (Sign Divergence Test)') ,"\n",
             '   Test selection format: ',ifelse(test == 'auto','Automatic selection','Manual selection') ,"\n","\n",
             '   Explanation:',"\n",
             '  -------------',"\n",
             ifelse(test == "auto",paste0(
               '   Automatic selection procedure based on the following:',"\n",
               ifelse(ttype == 'ENT',paste0('   The pooled effect (',round(rema_beta,4), ') was not statistically significant (p = ',round(rema_pval,4),')\n    or was exactly zero.'),
                      ifelse(ttype == 'SDT counter-positive',paste0('   The pooled effect (',round(rema_beta,4), ') was positive and statistically significant (p = ',round(rema_pval,4),').'),
                             paste0('   The pooled effect (',round(rema_beta,4), ') was negative and statistically significant (p = ',round(rema_pval,4),').')))),
               '   You manually selected the test.'),
             "\n",ifelse(ttype == 'ENT',paste0('   Testing whether number of all statistically significant effects (on both sides)\n    is greater than ',pnull),
                         ifelse(ttype == 'SDT counter-positive',paste0('   Testing whether number of negative statistically significant effects\n    is greater than ',pnull),
                                paste0('   Testing whether number of positive statistically significant effects\n    is greater than ',pnull))),"\n",
             "\n",'   Results:',"\n",
             '  ---------',"\n",
             '   Context - REMA pooled effect (p-val): ',round(rema_beta,3),'(',round(rema_pval,4),')',"\n","\n",
             '   Total number of significant effects (both sides): ',all_sig_sum,"\n",
             '   Number of significant positive cases: ',positive_sig_sum,"\n",
             '   Number of significant negative cases: ',negative_sig_sum,"\n",
             '   Number of valid cases: ',number_of_effects, "\n","\n",
             '   Given the test, relevant significant effects are ',ifelse(ttype == 'ENT', 'all (both sides): ',ifelse(ttype == 'SDT counter-positive', 'just negatives: ','just positives: ')), sig_effects,"\n",
             '   Number of valid cases: ',number_of_effects, "\n","\n",
             '   ',ifelse(ttype == 'ENT','ENT ','SDT '),'p-value = ', round(bintest$p.value,7)))


  #Create list of parameteres for return.

  sden_params <- list(test_type = ttype, selection_mechanism = test, rema_beta = rema_beta, rema_pval = rema_pval, all_sig_sum = all_sig_sum,
                      positive_sig_sum = positive_sig_sum, negative_sig_sum = negative_sig_sum, number_of_effects = number_of_effects, test_pval = bintest$p.value)

  return(list(sden_parameters = sden_params, binomial_test = bintest))
}

