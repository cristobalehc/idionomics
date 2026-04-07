#######################################
############ sden test     ############
#####################################

#' Run Sign Divergence or Equisyncratic Null tests.
#'
#' @export
#'
#' @param iarimax_object Your iarimax object.
#' @param alpha_arimax Significance threshold for classifying individual subjects'
#'   ARIMAX coefficients as significant (\code{pval < alpha_arimax}).
#'   Also used to derive the null probability for the binomial test:
#'   \code{alpha_arimax} in ENT mode, \code{alpha_arimax / 2} in SDT mode
#'   (unless overridden by \code{alpha_binom}). Defaults to 0.05.
#'   In ENT mode this creates a double role: raising \code{alpha_arimax}
#'   simultaneously increases the count of significant effects and raises the
#'   null expectation against which they are tested. Use \code{alpha_binom} to
#'   decouple the two.
#' @param alpha_binom Null probability for the binomial test. Defaults to \code{alpha_arimax} for
#'   ENT and \code{alpha_arimax / 2} for SDT. Supply this to override the default without changing
#'   the per-subject significance threshold.
#' @param test Type of test, default to "auto". "auto", "SDT", and "ENT" are supported.
#'   When `"auto"`, the selection between SDT and ENT is based on whether the pooled REMA
#'   effect is statistically significant at a **fixed threshold of 0.05**, regardless of
#'   `alpha_arimax`. This is intentional: the pooled-effect pivot and the per-subject
#'   significance criterion are treated as separate inferential decisions. Note that this
#'   threshold also appears in the advisory messages printed under `"SDT"` and `"ENT"`,
#'   but does not affect any computation in those modes.
#' @param feature Feature name to run sden test. Defaults to iarimax focal predictor. Use your original name, function will automatically append "estimate_"
#'
#' @return An S3 object of class \code{sden_results} with two elements:
#'   \code{sden_parameters} (a named list with test type, selection mechanism,
#'   REMA beta, p-null, significance counts, and binomial p-value) and
#'   \code{binomial_test} (the \code{htest} object from \code{stats::binom.test}).
#'   Attributes \code{focal_predictor}, \code{id_var}, and \code{timevar} are
#'   inherited from the input \code{iarimax_object}.
#'
#' @seealso [iarimax()] for per-subject modeling, [i_pval()] for per-subject
#'   p-values (called internally by `sden_test()`),
#'   [summary.sden_results()] for printing results.
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' panel <- do.call(rbind, lapply(1:6, function(id) {
#'   x <- rnorm(30)
#'   data.frame(id = as.character(id), time = seq_len(30),
#'              x = x, y = 0.5 * x + rnorm(30),
#'              stringsAsFactors = FALSE)
#' }))
#'
#' result <- iarimax(panel, y_series = "y", x_series = "x",
#'                   id_var = "id", timevar = "time",
#'                   min_n_subject = 20)
#'
#' # Automatic test selection based on pooled REMA effect
#' sden <- sden_test(result)
#' summary(sden)
#'
#' # Force ENT regardless of REMA
#' sden_ent <- sden_test(result, test = "ENT")
#' }

sden_test <- function(iarimax_object,alpha_arimax = 0.05, alpha_binom = NULL, test = "auto", feature = NULL){

  #Guard clause.

  #Only works for iarimax object now.
  if (!inherits(iarimax_object, "iarimax_results")) { #will add other classes if needed.
    stop("iarimax_object must be an iarimax_results object.")
  }

  # Check if test is valid
  if (!test %in% c('auto', 'SDT', 'ENT')) {
    stop('Invalid test type. Supported tests are "auto", "SDT", and "ENT".')
  }

  if (!is.numeric(alpha_arimax) || length(alpha_arimax) != 1 ||
      !is.finite(alpha_arimax) || alpha_arimax <= 0 || alpha_arimax >= 1) {
    stop("'alpha_arimax' must be a single finite number in (0, 1).")
  }
  if (!is.null(alpha_binom) &&
      (!is.numeric(alpha_binom) || length(alpha_binom) != 1 ||
       !is.finite(alpha_binom) || alpha_binom <= 0 || alpha_binom >= 1)) {
    stop("'alpha_binom' must be a single finite number in (0, 1), or NULL.")
  }

  #Defaults to focal predictor if feature is null.
  if (is.null(feature)) {
    feature <- attr(iarimax_object, "focal_predictor")
  }

  #Start the test.

  # Get p-values based on t-distribution.
  iarimax_object <- i_pval(iarimax_object, feature = feature)

  # Get number of valid effects (with p-values).
  pval_col_name <- paste0("pval_", feature) #Build pval column name.
  number_of_effects <- sum(!is.na(iarimax_object$results_df[[pval_col_name]]))

  #Create feature name.
  est_col_name <- paste0("estimate_", feature)
  std_col_name <- paste0("std.error_", feature)

  #Check if focal predictor is used, if not run RMA again.
  if (feature != attr(iarimax_object, "focal_predictor")) {


    #Override meta analysis object within the context of the function.
    iarimax_object$meta_analysis <-
      tryCatch(
        {
          metafor::rma(yi = iarimax_object$results_df[[est_col_name]], sei = iarimax_object$results_df[[std_col_name]], method = "REML")
        },
        error = function(e) {
          message('Error running meta-analysis: ', e$message, '. SDEN test will not be computed')
          NULL
        }
      )
  }

  #Stop if meta analysis is null.
  if (is.null(iarimax_object$meta_analysis)) {
    stop("sden_test requires a valid meta_analysis in the iarimax object, but it is NULL. ",
         "Check your iarimax output before calling sden_test().")
  }

  # Get REMA parameters.
  rema_beta <- as.numeric(iarimax_object$meta_analysis$beta)
  rema_pval <- iarimax_object$meta_analysis$pval

  # Get number of statistically signficant cases.
  positive_sig_sum <- sum(ifelse(iarimax_object$results_df[[est_col_name]] > 0 &
                                   iarimax_object$results_df[[pval_col_name]] < alpha_arimax,1,0), na.rm = TRUE)
  negative_sig_sum <- sum(ifelse(iarimax_object$results_df[[est_col_name]] < 0 &
                                   iarimax_object$results_df[[pval_col_name]] < alpha_arimax,1,0), na.rm = TRUE)
  all_sig_sum <- sum(ifelse(iarimax_object$results_df[[pval_col_name]] < alpha_arimax,1,0), na.rm = TRUE)


  #Automatically select the test.
  if (test == "auto") {

    if (rema_pval >= 0.05) {
      message('Running Equisyncratic Null Test (ENT)')

      sig_effects <-  all_sig_sum
      pnull <- if (!is.null(alpha_binom)) alpha_binom else alpha_arimax
      ttype <- 'ENT'

    }
    else if (rema_pval < 0.05 && rema_beta > 0) {
      message('Running Sign Divergence Test (SDT): number of negative cases (counter positive pooled-effect) values are being evaluated')

      sig_effects <- negative_sig_sum
      pnull <- if (!is.null(alpha_binom)) alpha_binom else alpha_arimax / 2
      ttype <- 'SDT counter-positive'

    }
    else if (rema_pval < 0.05 && rema_beta < 0) {

      message('Running Sign Divergence Test (SDT): number of positive cases (counter negative pooled-effect) values are being evaluated')

      sig_effects <- positive_sig_sum
      pnull <- if (!is.null(alpha_binom)) alpha_binom else alpha_arimax / 2
      ttype <- 'SDT counter-negative'

    }
    else {
      message('Meta-Analysis beta value is exactly zero, defaulting to ENT.')
      message('Running Equisyncratic Null Test (ENT).')

      sig_effects <- all_sig_sum
      pnull <- if (!is.null(alpha_binom)) alpha_binom else alpha_arimax
      ttype <- 'ENT'

    }

  }

  else if (test == "SDT") {

    if (rema_beta > 0) {
      message('Running Sign Divergence Test (SDT): number of negative cases (counter positive pooled-effect) values are being evaluated')

      if (rema_pval >= 0.05) {
        message('The p-value of the pooled effect is not significant at 0.05.',
                ' p = ', round(rema_pval, 4), '. ENT test is suggested.')
      }

      sig_effects <- negative_sig_sum
      pnull <- if (!is.null(alpha_binom)) alpha_binom else alpha_arimax / 2
      ttype <- 'SDT counter-positive'


    }
    else if (rema_beta < 0) {

      if (rema_pval >= 0.05) {
        message('The p-value of the pooled effect is not significant at 0.05.',
                ' p = ', round(rema_pval, 4), '. ENT test is suggested.')
      }

      message('Running Sign Divergence Test (SDT): number of positive cases (counter negative pooled-effect) values are being evaluated')

      sig_effects <- positive_sig_sum
      pnull <- if (!is.null(alpha_binom)) alpha_binom else alpha_arimax / 2
      ttype <- 'SDT counter-negative'


    }
    else {
      message('Meta-Analysis beta value is exactly zero, defaulting to ENT.')
      message('Running Equisyncratic Null Test (ENT).')

      sig_effects <- all_sig_sum
      pnull <- if (!is.null(alpha_binom)) alpha_binom else alpha_arimax
      ttype <- 'ENT'

    }

  }
  else if (test == "ENT") {

    message('Running Equisyncratic Null Test (ENT)')

    if (rema_pval < 0.05) {
      message('The p-value of the pooled effect is statistically significant at 0.05.',
              ' p = ', round(rema_pval, 4), '. SDT test is suggested.')
    }

    sig_effects <- all_sig_sum
    pnull <- if (!is.null(alpha_binom)) alpha_binom else alpha_arimax
    ttype <- 'ENT'

  }
  #Run the appropriate binomial test.

  #Guard: binom.test requires n > 0.
  if (number_of_effects == 0) {
    stop('No valid p-values found for feature "', feature, '". All subjects may have failed model fitting.')
  }

  bintest <- stats::binom.test(x = sig_effects, n = number_of_effects, p = pnull, alternative = "greater")


  #Create list of parameters for return.

  sden_params <- list(test_type = ttype, selection_mechanism = test, rema_beta = rema_beta, pnull = pnull, rema_pval = rema_pval, all_sig_sum = all_sig_sum,
                      positive_sig_sum = positive_sig_sum, negative_sig_sum = negative_sig_sum, number_of_effects = number_of_effects, test_pval = bintest$p.value,
                      sig_effects = sig_effects)


  result <- list(sden_parameters = sden_params, binomial_test = bintest)

  #Add class for S3 dispatch.
  class(result) <- c("sden_results","list")

  #Inherit attributes to identify focal predictor, and id_var.
  attr(result, "outcome") <- attr(iarimax_object, "outcome")
  attr(result, "focal_predictor") <- feature
  attr(result, "id_var") <- attr(iarimax_object, "id_var")
  attr(result, "timevar") <- attr(iarimax_object, "timevar")


  #Return.
  return(invisible(result))
}
