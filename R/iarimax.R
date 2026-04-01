#' Run I-ARIMAX algorithm.
#'
#' @export
#' @importFrom utils globalVariables
#' @importFrom rlang :=
#'
#' @param dataframe Your dataframe.
#' @param min_n_subject Integer. Subjects with fewer than \code{min_n_subject}
#'   pairwise-complete observations (across \code{y_series} and all
#'   \code{x_series}) are excluded. The threshold is inclusive (\code{>=}).
#'   \code{n_valid} is extracted directly from the fitted Arima object via
#'   \code{stats::nobs()} and might be smaller than \code{min_n_subject} if
#'   \code{auto.arima()} applies differencing. Defaults to 20.
#' @param minvar Numeric. Each series (\code{y_series} and every variable in
#'   \code{x_series}) must have variance \code{>= minvar}; subjects failing
#'   this for any series are excluded. Defaults to 0.01. Person-mean standardizing
#'   variables will make this parameter irrelevant, as var == 1 under z standardization.
#' @param y_series A string containing the name of your dependent variable y.
#' @param x_series A character vector containing the name(s) of your predictor variable(s).
#' @param focal_predictor A string with the name of the predictor to use in the meta-analysis.
#'  Required when x_series has more than one variable. When x_series is a single variable, defaults to that variable.
#' @param id_var A string containing your id variable.
#' @param timevar A string naming the column used to sort observations into
#'   chronological order within each subject. Must be complete (no missing
#'   values) and numeric.
#' @param fixed_d Optional non-negative integer. If provided, fixes the
#'   differencing order to this value for every subject instead of letting
#'   \code{auto.arima()} select it. \code{NULL} (default) means automatic
#'   selection. Fixing \code{d} ensures all subjects' coefficients are on
#'   the same scale (e.g., \code{d = 0} for levels, \code{d = 1} for
#'   changes), which is important for comparability in the meta-analysis.
#'   AR (\code{p}) and MA (\code{q}) orders are always selected
#'   automatically per subject because they control error autocorrelation
#'   and do not affect the scale of the xreg coefficients.
#' @param correlation_method Select method for raw correlations. Options are: 'spearman', 'pearson' or 'kendall'. Defaults to 'pearson'.
#' @param keep_models If TRUE, will keep original arimax models in a list.
#' @param verbose If TRUE, prints progress messages during filtering and model fitting. Defaults to FALSE.
#'
#' @return An S3 object of class \code{iarimax_results} (a named list) with:
#' \describe{
#'   \item{results_df}{Per-subject data frame with ARIMA orders (\code{nAR},
#'     \code{nI}, \code{nMA}), coefficient estimates
#'     (\code{estimate_<feature>}), standard errors
#'     (\code{std.error_<feature>}), \code{n_valid} (effective n extracted
#'     from the fitted Arima object via \code{stats::nobs()}),
#'     \code{n_params} (number of model coefficients), and \code{raw_cor}
#'     (pairwise correlation between the focal predictor and
#'     \code{y_series}, not partialled for other predictors in
#'     \code{x_series}). Subjects whose \code{auto.arima()} call failed
#'     appear with \code{NA} for all model statistics (\code{nAR},
#'     \code{nI}, \code{nMA}, \code{n_valid}, \code{n_params}, and all
#'     coefficient/SE columns), but retain a valid \code{raw_cor} if the
#'     correlation step succeeded before the model failure. \code{raw_cor}
#'     is \code{NA} only if \code{cor.test()} also failed. Their IDs are
#'     listed in \code{case_number_detail$error_arimax_skipped} if auto.arima fails.}
#'   \item{meta_analysis}{A \code{metafor::rma} object from a random-effects
#'     meta-analysis on the focal predictor coefficients, or \code{NULL} if
#'     the meta-analysis failed (e.g. fewer than two valid estimates).}
#'   \item{case_number_detail}{A list with \code{n_original_df} (total unique
#'     subjects in the input), \code{n_filtered_out} (excluded by variance or
#'     n thresholds), \code{error_arimax_skipped} (character vector of IDs
#'     whose \code{auto.arima()} call failed), and \code{n_used_iarimax}
#'     (subjects with a valid fitted model). The four quantities satisfy:
#'     \code{n_original_df == n_filtered_out + length(error_arimax_skipped) +
#'     n_used_iarimax}.}
#'   \item{models}{A named list of raw \code{Arima} objects (one per subject,
#'     \code{NULL} for failed fits) if \code{keep_models = TRUE}, otherwise
#'     \code{NULL}.}
#' }
#' Attributes: \code{outcome} (the \code{y_series} name),
#' \code{focal_predictor}, \code{id_var}, \code{timevar}.
#'
#' @examples
#' \donttest{
#' # Build a small synthetic panel: 6 subjects, 30 observations each
#' set.seed(42)
#' panel <- do.call(rbind, lapply(1:6, function(id) {
#'   x <- rnorm(30)
#'   data.frame(
#'     id   = as.character(id),
#'     time = seq_len(30),
#'     x    = x,
#'     y    = 0.5 * x + rnorm(30),
#'     stringsAsFactors = FALSE
#'   )
#' }))
#'
#' result <- iarimax(panel,
#'                   y_series  = "y",
#'                   x_series  = "x",
#'                   id_var    = "id",
#'                   timevar   = "time")
#'
#' summary(result)
#' plot(result)
#' }


#######################################
############ I ARIMAX FUNCTION #######
#####################################

iarimax <- function(dataframe, min_n_subject = 20, minvar = 0.01, y_series, x_series,
                    focal_predictor = NULL, id_var, timevar, fixed_d = NULL,
                    correlation_method = 'pearson', keep_models = FALSE, verbose = FALSE) {

  # Check whether variables are in the dataset.
  required_vars <- c(y_series, x_series, id_var, timevar)

  if (!is.numeric(min_n_subject) || length(min_n_subject) != 1 ||
      !is.finite(min_n_subject) || min_n_subject < 1) {
    stop("'min_n_subject' must be a finite positive number.")
  }
  if (!is.numeric(minvar) || length(minvar) != 1 ||
      !is.finite(minvar) || minvar < 0) {
    stop("'minvar' must be a finite non-negative number.")
  }

  if (!is.null(fixed_d)) {
    if (!is.numeric(fixed_d) || length(fixed_d) != 1 ||
        !is.finite(fixed_d) || fixed_d < 0 || fixed_d != round(fixed_d)) {
      stop("'fixed_d' must be a single non-negative integer (e.g., 0 or 1).")
    }
    fixed_d <- as.integer(fixed_d)
  }

  if (!correlation_method %in% c("pearson","spearman","kendall")) {
    stop(paste("Correlation method not supported. Check if you spelled the method correctly:", correlation_method))
  }

  if (!all(required_vars %in% colnames(dataframe))) {
    missing_vars <- required_vars[!required_vars %in% colnames(dataframe)]
    stop(paste("Cannot find required variables. Check if you spelled the following variables correctly:", paste(missing_vars, collapse = ", ")))
  }

  # timevar must be numeric for correct temporal ordering
  if (!is.numeric(dataframe[[timevar]])) {
    stop(
      "Column '", timevar, "' must be numeric. Got class: ",
      class(dataframe[[timevar]])[1], ". Convert dates or other formats to a reasonable numeric value before calling iarimax."
    )
  }

  #Check for timevar missing data: With missings, the data cannot be arranged.
  n_missing_timevar <- sum(is.na(dataframe[[timevar]]))
  if (n_missing_timevar > 0) {
    stop(
      n_missing_timevar, " row(s) have missing values in the time variable '", timevar, "'. ",
      "iarimax requires a non-missing timevar for every observation to ensure correct temporal ordering. ",
      "Please review these rows and process accordingly before running iarimax."
    )
  }

  # Resolve focal_predictor: default to x_series when only one predictor is given.
  if (!is.null(focal_predictor)) {
    if (!is.character(focal_predictor) || length(focal_predictor) != 1) {
      stop("'focal_predictor' must be a single character string naming one of the x_series variables.")
    }
  }

  if (is.null(focal_predictor)) {
    if (length(x_series) == 1) {
      focal_predictor <- x_series
    } else {
      stop("focal_predictor is required when x_series contains more than one variable.")
    }
  }

  if (!focal_predictor %in% x_series) {
    stop(paste("focal_predictor must be one of the x_series variables. Got:", focal_predictor))
  }

  #Extract number of id's.
  number_of_original_ids <- length(unique(dataframe[[id_var]]))

  # Convert strings to rlang::symbols
  y_series_sym <- rlang::sym(y_series)
  id_var_sym <- rlang::sym(id_var)
  timevar_sym <- rlang::sym(timevar)

  #Id var as character, to avoid numerics and giant lists.
  dataframe[[id_var]] <- as.character(dataframe[[id_var]])


  #Subset only relevant features.
  dataframe <- dataframe[, required_vars, drop = FALSE]

  if (verbose) message('Filtering data based on minimum non-NA observations and variance...')

  # Filter N listwise complete observations with variance conditions
  subjects <- dataframe |>
    dplyr::group_by(!!id_var_sym) |> #Group by id variable.
    dplyr::filter(!is.na(!!y_series_sym) &
                    dplyr::if_all(dplyr::all_of(x_series), ~ !is.na(.x))) |> #Filter all that are complete in all variables.
    dplyr::summarise(
      count = dplyr::n(), #Use counts to filter.
      var_outcome = stats::var(!!y_series_sym, na.rm = TRUE), #use variance to filter.
      dplyr::across(dplyr::all_of(x_series), ~ stats::var(.x, na.rm = TRUE), .names = "var_{.col}"), #Variance for each predictor.
      .groups = 'drop'
    ) |>
    dplyr::filter(count >= min_n_subject,
                  dplyr::if_all(dplyr::starts_with("var_"), ~ .x >= minvar)) |> #y and all predictors must meet minvar.
    dplyr::pull(!!id_var_sym) #Filter data.

  #ID as character, to avoid giant lists.
  subjects <- as.character(subjects)

  #Stop if not enough subjects.
  if (length(subjects) <= 1){
    stop(paste("There are not enough cases to run the iarimax algorithm. You have", length(subjects),"valid cases."))
  }

  if (verbose) {
    message('Filtering done: ', length(subjects), ' subjects will be used for the analyses.')
    message('Running I-ARIMAX algorithm...')
  }

  #Storage lists.

  #Number of valid cases & parameters.
  n_valid <- list()
  n_params <- list()

  #Exclude cases where arimax don't work.
  exclude <- character(0)

  #Model lists.
  arimax_models <- list()
  tidy_models <- list()
  raw_correlation <- list()
  AR_N <- list()
  I_N <- list()
  MA_N <- list()


  #Start case number counter.
  casen <- 0
  for(i in subjects) {

    #Update case number.
    casen <- casen + 1

    if (verbose) {
      d_label <- if (is.null(fixed_d)) "auto" else paste0("d=", fixed_d)
      message('  Applying ARIMAX (', d_label, ') to case: ', i, ' ... ', appendLF = FALSE)
    }

    #Extract the current subject & arrange timeseries by timevar.
    subject_n <- dataframe |>
      dplyr::filter(!!id_var_sym == i) |>
      dplyr::arrange(!!timevar_sym) #Ensure time-series order.

    ##########################################
    ### A COMMENT OF MISSING DATA HANDLING ###
    #########################################

    # Note: We do NOT need to manually sync NAs.
    # stats::arima's Kalman Filter (C source) automatically treats
    # any row with a missing predictor as a missing observation
    # in the state-space update, preserving temporal structure.

    #Extract y vector.
    y_vector <- subject_n |>
      dplyr::pull(!!y_series_sym)

    #Extract x matrix (one column per predictor, named).
    x_matrix <- as.matrix(subject_n[, x_series, drop = FALSE])


    ########################
    #### Calculate cor #####
    ########################

    # Correlation is computed for the focal predictor only.
    focal_vector <- subject_n |> dplyr::pull(focal_predictor)

    correlation <- tryCatch(
      {#Supress spearman's warning about p-values with ties.
        suppressWarnings(stats::cor.test(x = focal_vector, y = y_vector, method = correlation_method))
      },
      error = function(e) {
        message('\nError computing correlation for case: ', i, '\n  ', e$message)
        NULL #Set model as NULL
      }
    )

    #Handle when correlations are null.
    if (is.null(correlation)) {
      raw_correlation[[i]] <- NA
    } else {
      raw_correlation[[i]] <- correlation$estimate[1]
    }

    ##############################
    ### CALCULATE AUTO.ARIMA ####
    ############################

    model <- tryCatch(
      {
        forecast::auto.arima(y = y_vector, xreg = x_matrix, d = if (is.null(fixed_d)) NA else fixed_d, approximation = FALSE, stepwise = FALSE)
      },
      error = function(e) {
        message('\nError running ARIMAX model for case: ', i, '\n  ', e$message)
        NULL #Set model as NULL
      }
    )

    #Fill the list with NA if the model is null

    if (is.null(model)) {
      AR_N[[i]] <- NA
      I_N[[i]] <- NA
      MA_N[[i]] <- NA
      n_valid[[i]] <- NA
      n_params[[i]] <- NA
      tidy_models[[i]] <- NULL
      arimax_models[[i]] <- NULL

      exclude <- c(exclude,i)

      if (verbose) message('skipped. (', round(casen / length(subjects) * 100, 1), '% done)')
      next #Skip the next part of this iteration of the loop, so it doesn't get overriden and throws an error.
    }

    # Remove fable's "ARIMA" class to prevent it from hijacking broom::tidy.
    class(model) <- setdiff(class(model), "ARIMA")
    tidymodel <- broom::tidy(model)

    #Add id to tidy model.
    tidymodel <- tidymodel |>
      dplyr::mutate(!!id_var_sym := i)

    #Fill the tidy model list.
    tidy_models[[i]] <- tidymodel
    # Conditional storage of raw models
    if (keep_models) {
      arimax_models[[i]] <- model
    } else {
      arimax_models[[i]] <- NULL # Placeholder to keep list lengths consistent
    }


    #Fill the number of ARIMA parameteres lists: Just the number of AR I MA processes involved.
    AR_N[[i]] <- model$arma[1] #Fill AR list.
    I_N[[i]] <- model$arma[6] #Fill I list.
    MA_N[[i]] <- model$arma[2] #Fill MA List.

    #Fill n valid and n params. nobs() returns the n used in the model likelihood.
    n_valid[[i]] <- stats::nobs(model)
    n_params[[i]] <- length(model$coef)


    if (verbose) message(round(casen / length(subjects) * 100, 1), '% done')
  }

  ###############################################
  ### CREATE DATAFRAME WITH PARAMETERS ##########
  ###############################################
  #Filter null models.
  tidy_list <- Filter(Negate(is.null), tidy_models)

  #Stop if no valid model is present.
  if (length(tidy_list) == 0) {
    stop(
      "iarimax: all ARIMAX fits failed after filtering. ",
      "No subject produced a valid model. ",
      "Check (a) time ordering / timevar, (b) missingness in y/x, ",
      "(c) min_n_subject/minvar thresholds, and (d) whether xreg has NA patterns."
    )
  }

  #Create tidy long.
  tidy_long <- dplyr::bind_rows(tidy_list)


  # Pivot the coefficients to wide format.
  # Term names come from the xreg column names, e.g. estimate_stress, std.error_stress, estimate_ar1, etc.
  tidy_wide <- tidy_long |>
    tidyr::pivot_wider(
      id_cols     = !!id_var_sym,
      names_from  = "term",
      values_from = c("estimate", "std.error")
    )


  # Convert the lists to vectors
  AR_vector <- unlist(AR_N)
  I_vector <- unlist(I_N)
  MA_vector <- unlist(MA_N)
  n_valid_vector <- unlist(n_valid)
  n_params_vector <- unlist(n_params)
  raw_correlation_vector <- unlist(raw_correlation)


  #Create individual level summaries.
  summary_df <- tibble::tibble(
    !!id_var_sym := subjects,
    nAR = AR_vector,
    nI = I_vector,
    nMA = MA_vector,
    raw_cor = raw_correlation_vector,
    n_valid = n_valid_vector,
    n_params = n_params_vector
  )

  #Create dataframe to return.

  results_df <- summary_df |>
    dplyr::left_join(tidy_wide, by = id_var)


  ##############################################
  ##### RUN RANDOM EFFECTS META ANALYSIS ######
  #############################################

  #Try running the random effects meta analysis.

  focal_est_col <- paste0("estimate_",  focal_predictor)
  focal_se_col  <- paste0("std.error_", focal_predictor)

  meta_analysis <-
    tryCatch(
      {
        metafor::rma(yi = results_df[[focal_est_col]], sei = results_df[[focal_se_col]], method = "REML")
      },
      error = function(e) {
        message('Error running meta-analysis: ', e$message, '. meta_analysis will be NULL in the output.')
        NULL
      }
    )

  final_obj <- list(results_df = results_df,
                    meta_analysis = meta_analysis,
                    case_number_detail = list(n_original_df = number_of_original_ids, #Number of original ids.
                                                n_filtered_out = number_of_original_ids-length(subjects), #Filtered out by var or n.
                                                error_arimax_skipped = exclude, #auto.arima failed.
                                                n_used_iarimax = nrow(results_df)-length(exclude)), #Final number of used cases.
                    models = if (keep_models) arimax_models else NULL)

  #Add class to identify the object.
  class(final_obj) <- c("iarimax_results", "list")

  #Add attribute to identify focal predictor, and id_var.
  attr(final_obj, "outcome") <- y_series
  attr(final_obj, "focal_predictor") <- focal_predictor
  attr(final_obj, "id_var") <- id_var
  attr(final_obj, "timevar") <- timevar


  if (verbose) message('I-ARIMAX algorithm finished.')

  return(final_obj)

}


utils::globalVariables(c("count", "var_outcome")) #Declare symbolic global variables.
