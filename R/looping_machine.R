#' Run dynamic looping I-ARIMAX algorithm.
#'
#' Fits three I-ARIMAX models forming a directed loop (a -> b -> c -> a),
#' applies per-subject p-values, and flags subjects where all three focal
#' coefficients are positive and significant.
#'
#' @export
#'
#' @param dataframe Your dataframe.
#' @param a_series,b_series,c_series Strings naming the three loop variables.
#' @param id_var String naming the subject ID variable.
#' @param timevar String naming the time variable (must be complete, no NAs).
#' @param covariates Optional character vector of additional predictors added to all three legs. Defaults to NULL.
#' @param include_third_as_covariate If TRUE, the third loop variable is added as a covariate in each leg. Note: this also applies \code{minvar} filtering to the third variable, which may silently reduce the subject count. Defaults to FALSE.
#' @param min_n_subject Minimum pairwise-complete observations per subject. Defaults to 20.
#' @param minvar Minimum variance required in all series to retain a subject. Defaults to 0.01.
#' @param correlation_method Raw correlation method: 'pearson', 'spearman', or 'kendall'. Defaults to 'pearson'.
#' @param alpha Significance threshold for the loop condition; must be in (0, 1). Defaults to 0.05.
#' @param keep_models If TRUE, raw ARIMAX objects are kept in each leg's result. Defaults to FALSE.
#' @param verbose If TRUE, prints per-leg progress messages. Defaults to FALSE.
#'
#' @section Statistical notes:
#'
#' **Conjunction criterion and Type I error.** \code{Loop_positive_directed} requires
#' all three tests to be significant at \code{alpha} and all three coefficients to
#' be positive. By Boole's inequality,
#' \eqn{P(T_1 \cap T_2 \cap T_3) \leq \min(P(T_i)) = \alpha} under the global
#' null, regardless of the correlation between legs. The criterion is always
#' conservative or at most exactly \eqn{\alpha}.
#'
#' **Partial effects on a possibly differenced outcome.** Coefficients are partial
#' regression effects conditioned on the ARMA structure chosen by
#' \code{auto.arima()}. When differencing is selected (\eqn{d > 0}), the
#' coefficient describes the relationship on the differenced scale, not the raw
#' level. Check the ARIMA orders in each leg's \code{results_df} before
#' interpreting the sign of \code{Loop_positive_directed}.
#'
#' @return A named list:
#' \describe{
#'   \item{loop_df}{Per-subject dataframe with coefficients, SEs, n_valid, n_params, p-values, and \code{Loop_positive_directed} (1 = positive loop, 0 = not, NA = ARIMAX failed in at least one leg).}
#'   \item{alpha}{Significance threshold used.}
#'   \item{covariates}{Covariate vector (NULL if none).}
#'   \item{include_third_as_covariate}{Logical; whether the third variable was added as a covariate.}
#'   \item{loop_case_detail}{List with \code{n_in_loop_df}, \code{n_complete}, and \code{n_na_indicator} (subjects with a failed model in any leg).}
#'   \item{iarimax_a_to_b,iarimax_b_to_c,iarimax_c_to_a}{\code{iarimax_results} objects for each leg, with \code{i_pval} already applied.}
#' }
#'
#' @examples
#' \donttest{
#' # Build a panel with three correlated variables
#' set.seed(7)
#' panel <- do.call(rbind, lapply(1:6, function(id) {
#'   a <- rnorm(30)
#'   b <- 0.4 * a + rnorm(30)
#'   c <- 0.4 * b + rnorm(30)
#'   data.frame(id = as.character(id), time = seq_len(30),
#'              a = a, b = b, c = c, stringsAsFactors = FALSE)
#' }))
#'
#' loop_result <- looping_machine(panel,
#'                                a_series = "a", b_series = "b", c_series = "c",
#'                                id_var   = "id", timevar  = "time")
#'
#' # Proportion of subjects with a detected positive directed loop
#' mean(loop_result$loop_df$Loop_positive_directed, na.rm = TRUE)
#' }


##############################################
############ Looping Machine function ########
##############################################

looping_machine <- function(dataframe, a_series, b_series, c_series, id_var, timevar,
                            covariates = NULL,
                            include_third_as_covariate = FALSE,
                            min_n_subject = 20, minvar = 0.01,
                            correlation_method = 'pearson',
                            alpha = 0.05,
                            keep_models = FALSE,
                            verbose = FALSE) {

  # Guard: all three series must be different variables.
  series_names <- c(a_series, b_series, c_series)
  if (length(unique(series_names)) < 3) {
    stop("a_series, b_series, and c_series must all refer to different variables. Got: ",
         paste(series_names, collapse = ", "))
  }

  # Guard: alpha must be a single finite value strictly in (0, 1).
  if (!is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single finite numeric value strictly between 0 and 1. Got: ", alpha)
  }

  # Guard: covariates must not overlap with the loop variables.
  if (!is.null(covariates)) {
    overlap <- intersect(covariates, series_names)
    if (length(overlap) > 0) {
      stop("'covariates' must not overlap with a_series, b_series, or c_series. ",
           "Overlapping variables: ", paste(overlap, collapse = ", "))
    }
  }

  # Create name tags.
  ab_name <- paste0(a_series, "_", b_series)
  bc_name <- paste0(b_series, "_", c_series)
  ca_name <- paste0(c_series, "_", a_series)

  # Build x_series for each leg, optionally including the third loop variable as covariate.
  x_ab <- c(a_series, covariates, if (include_third_as_covariate) c_series)
  x_bc <- c(b_series, covariates, if (include_third_as_covariate) a_series)
  x_ca <- c(c_series, covariates, if (include_third_as_covariate) b_series)

  # a to b leg.
  if (verbose) message('Calculating a to b: From ', a_series, ' to ', b_series, '...')
  a_to_b <- iarimax(dataframe, min_n_subject = min_n_subject, minvar = minvar,
                    y_series = b_series, x_series = x_ab, focal_predictor = a_series,
                    id_var = id_var, timevar = timevar,
                    correlation_method = correlation_method,
                    keep_models = keep_models, verbose = verbose)
  a_to_b <- i_pval(a_to_b, feature = a_series)

  # b to c leg.
  if (verbose) message('Calculating b to c: From ', b_series, ' to ', c_series, '...')
  b_to_c <- iarimax(dataframe, min_n_subject = min_n_subject, minvar = minvar,
                    y_series = c_series, x_series = x_bc, focal_predictor = b_series,
                    id_var = id_var, timevar = timevar,
                    correlation_method = correlation_method,
                    keep_models = keep_models, verbose = verbose)
  b_to_c <- i_pval(b_to_c, feature = b_series)

  # c to a leg.
  if (verbose) message('Calculating c to a: From ', c_series, ' to ', a_series, '...')
  c_to_a <- iarimax(dataframe, min_n_subject = min_n_subject, minvar = minvar,
                    y_series = a_series, x_series = x_ca, focal_predictor = c_series,
                    id_var = id_var, timevar = timevar,
                    correlation_method = correlation_method,
                    keep_models = keep_models, verbose = verbose)
  c_to_a <- i_pval(c_to_a, feature = c_series)

  # Extract focal predictor columns by name for each leg.
  cols_ab <- c(id_var, paste0("estimate_", a_series), paste0("std.error_", a_series),
               "n_valid", "n_params", paste0("pval_", a_series))
  a_to_b_sub <- a_to_b$results_df[, cols_ab]
  colnames(a_to_b_sub) <- c(id_var, ab_name, paste0("stderr_", ab_name),
                             paste0(ab_name, "_n_valid"), paste0(ab_name, "_n_params"),
                             paste0(ab_name, "_pval"))

  cols_bc <- c(id_var, paste0("estimate_", b_series), paste0("std.error_", b_series),
               "n_valid", "n_params", paste0("pval_", b_series))
  b_to_c_sub <- b_to_c$results_df[, cols_bc]
  colnames(b_to_c_sub) <- c(id_var, bc_name, paste0("stderr_", bc_name),
                             paste0(bc_name, "_n_valid"), paste0(bc_name, "_n_params"),
                             paste0(bc_name, "_pval"))

  cols_ca <- c(id_var, paste0("estimate_", c_series), paste0("std.error_", c_series),
               "n_valid", "n_params", paste0("pval_", c_series))
  c_to_a_sub <- c_to_a$results_df[, cols_ca]
  colnames(c_to_a_sub) <- c(id_var, ca_name, paste0("stderr_", ca_name),
                             paste0(ca_name, "_n_valid"), paste0(ca_name, "_n_params"),
                             paste0(ca_name, "_pval"))

  # Merge all three legs.
  loop_df <- merge(a_to_b_sub, b_to_c_sub, by = id_var)
  loop_df <- merge(loop_df, c_to_a_sub, by = id_var)

  # Create the positive directed loop indicator.
  # Subjects where any leg's ARIMAX failed will have NA pvals/estimates, so
  # Loop_positive_directed is NA for those subjects — a failed model is neither
  # confirmed as a loop nor confirmed as a non-loop.
  loop_df[["Loop_positive_directed"]] <- ifelse(
    loop_df[[paste0(ab_name, "_pval")]] < alpha &
      loop_df[[paste0(bc_name, "_pval")]] < alpha &
      loop_df[[paste0(ca_name, "_pval")]] < alpha &
      loop_df[[ab_name]] > 0 &
      loop_df[[bc_name]] > 0 &
      loop_df[[ca_name]] > 0,
    1L, 0L
  )

  if (verbose) message('Looping Machine finished.')

  n_na_indicator <- sum(is.na(loop_df[["Loop_positive_directed"]]))
  n_complete     <- sum(!is.na(loop_df[["Loop_positive_directed"]]))

  message(
    'Number of cases with the positive loop present: ',
    sum(loop_df[["Loop_positive_directed"]], na.rm = TRUE),
    if (n_na_indicator > 0)
      paste0(' (', n_na_indicator,
             ' subject(s) excluded: incomplete model fit in at least one leg)')
    else ''
  )

  return(list(
    loop_df                    = loop_df,
    alpha                      = alpha,
    covariates                 = covariates,
    include_third_as_covariate = include_third_as_covariate,
    loop_case_detail           = list(
      n_in_loop_df   = nrow(loop_df),
      n_complete     = n_complete,
      n_na_indicator = n_na_indicator
    ),
    iarimax_a_to_b             = a_to_b,
    iarimax_b_to_c             = b_to_c,
    iarimax_c_to_a             = c_to_a
  ))

}
