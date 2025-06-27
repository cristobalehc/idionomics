#' Create idionomic plots.
#'
#' @export
#' @importFrom utils globalVariables
#'
#' @param iarimax_object Your iarimax object.
#' @param id_var Your id variable.
#' @param y_series_name Optional: A string containing the substantive name of your y series to be used in the plot title.
#' @param x_series_name Optional: A string containing the substantive name of your x series to be used in the plot title.
#' @param alpha_crit_t Critical value for model significance. Defaluts to 0.05.
#' @param plot_type Plot type. Supports 'caterpillar','ecdf', and 'density'. Defaults to 'caterpillar'.
#' @param lims a vector of length 2 with limits for the density plot. Defaults to c(-1,1)
#'
#' @returns The plots object.


#######################################
############ PLOTOMATRON #######
#####################################



Plotomatron <- function(iarimax_object,id_var ,y_series_name = NULL, x_series_name = NULL, alpha_crit_t = 0.05, plot_type = 'caterpillar', lims = c(-1,1)) {


  #Guard clause.

  # Check if plot_type is valid
  if (!plot_type %in% c('caterpillar', 'ecdf', 'density')) {
    stop('Invalid plot_type. Supported plot types are "caterpillar", "ecdf", and "density".')
  }


  # Check if arimax_object is an iarimax object
  if (is.list(iarimax_object) && all(c('results_df', 'meta_analysis') %in% names(iarimax_object))) {
    df_plt <- iarimax_object$results_df
  } else {
    stop('Plotomatron requieres an iarimaxoid object that contains "results_df" and "meta_analysis". Run IARIMAXoid_Pro with metaanalysis = TRUE')
  }


  #Plot for contemporaneous relationships:


    #Create symbolic variables.
    id_var_sym <- rlang::sym(id_var)


  #Calcualte degrees of freedom and critical value.
  df_plt <- df_plt %>%
    dplyr::mutate(df_mod = n_valid - n_params)
  #Critical value.
  crit_val <- stats::qt(1 - (alpha_crit_t / 2), df_plt$df_mod)


  ### PLOTS FOR CONTEMPORANEOUS RELATIONSHIPS.

  if(iarimax_object$type == 'contemporaneous_predictor'){

  #Create Caterpillar Plot.
  if (plot_type == 'caterpillar') {
    # Reorder data
    df_plt <- df_plt %>%
      dplyr::mutate(!!id_var_sym := forcats::fct_reorder(as.factor(!!id_var_sym), xreg))

    # Create color coding
    df_plt <- df_plt %>%
      dplyr::mutate(line_color = dplyr::case_when(
        xreg - crit_val * stderr_xreg > 0 ~ "green",  # Interval entirely positive
        xreg + crit_val * stderr_xreg < 0 ~ "red",    # Interval entirely negative
        TRUE ~ "black"                                # Interval crosses zero or includes it
      ))

    # Create the plot
    caterpillar_plot <- ggplot2::ggplot(df_plt, ggplot2::aes(x = !!id_var_sym, y = xreg)) +
      ggplot2::geom_point() +  # Plot the mean intercepts as points
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::geom_linerange(ggplot2::aes(ymin = xreg - crit_val * stderr_xreg, ymax = xreg + crit_val * stderr_xreg, color = line_color)) +  # Add error bars for +/- 2SD with color based on the new variable
      ggplot2::scale_color_identity() + # Use the colors specified in the iarimax_object
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +  # Add a horizontal line at y = 0 for reference
      ggplot2::coord_flip() +  # Flip the coordinates to make it a horizontal plot
      ggplot2::labs(x = "Participant ID",
                    y = "i-ARIMAX Effect Sizes and 95% Confidence Intervals of \nWithin-Person Associations",
                    title = ifelse(is.null(y_series_name) | is.null(x_series_name), 'X series Linked to Y series', paste0(x_series_name, " Linked to  \n", y_series_name)),
                    caption = "Note: Blue vertical lines represent the RE-MA pooled (nomothetic) effect in the middle\n and the lower and upper bounds of the 95% CI of the pooled effect.\nGreen horizontal lines indicate 95% CI of positive associations and red indicate negative.") +
      ggplot2::geom_hline(yintercept = iarimax_object$meta_analysis$beta, color = "blue", linewidth = 1) +  # Add population-level mean intercept line
      ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = iarimax_object$meta_analysis$ci.lb, ymax = iarimax_object$meta_analysis$ci.ub, alpha = 0.2, color = "blue") + # Add 95% credible interval
      ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0)) +
      ggplot2::theme(axis.text.y = ggplot2::element_blank())

    return(caterpillar_plot)

  }

  #Plot empirical cumulative distribution function.
  if (plot_type == "ecdf") {

    if('rand_df' %in% names(iarimax_object)) {
      #Merge random effects to I-arimax dataframe.
      df_plt <- merge(df_plt, subset(iarimax_object$rand_df, select = c(id_var,'random_slope')), by = id_var, all.x=TRUE)
    }
    else {
      stop('Random effects are requiered for ecdf plot. Please run IARIMAXoid_Pro with hlm_compare = TRUE.')
      }

    #Pivot dataset to long format.
    df_plt_long <- df_plt %>%
      dplyr::select(c(!!id_var_sym,"xreg","random_slope","raw_correlation")) %>%
      dplyr::rename(Iarimax = xreg, MLM = random_slope, RawCorr = raw_correlation) %>%
      tidyr::pivot_longer(
        cols = c("Iarimax","MLM","RawCorr"),
        names_to = "Method",
        values_to = "Estimate"
      )

    #Create plot.
    ecdf_plot <- ggplot2::ggplot(df_plt_long, ggplot2::aes(Estimate, colour = Method)) +
      ggplot2::stat_ecdf(geom = "step", pad = FALSE, linewidth = 1) +
      ggplot2::ylab("Empirical Cumulative Distribution Function") +
      ggplot2::xlab(ifelse(is.null(y_series_name) | is.null(x_series_name),"Effect Size of the Link Between \n X series and Y series in Daily Life",
                           paste0("Effect Size of the Link Between \n",x_series_name, " and ",y_series_name, " in Daily Life"))) +
      ggplot2::scale_y_continuous(breaks = seq(0.00, 1.00, by = 0.10)) +
      ggplot2::scale_x_continuous(breaks = seq(-1.0, 1.0, by = 0.20)) +
      ggplot2::theme_bw()

    return(ecdf_plot)

  }

  if (plot_type == "density") {

    if('rand_df' %in% names(iarimax_object)) {
      #Merge random effects to I-arimax dataframe.
      df_plt <- merge(df_plt, subset(iarimax_object$rand_df, select = c(id_var,'random_slope')), by = id_var, all.x=TRUE)
    }
    else {
      stop('Random effects are requiered for density plot. Please run IARIMAXoid_Pro with hlm_compare = TRUE.')
    }

    #Pivot dataset to long format.
    df_plt_long <- df_plt %>%
      dplyr::select(c(!!id_var_sym,"xreg","random_slope","raw_correlation")) %>%
      dplyr::rename(Iarimax = xreg, MLM = random_slope, RawCorr = raw_correlation) %>%
      tidyr::pivot_longer(
        cols = c("Iarimax","MLM","RawCorr"),
        names_to = "Method",
        values_to = "Estimate"
      )

    #Plot density plot.
    density_plot <- df_plt_long %>%
      ggplot2::ggplot(ggplot2::aes(x = Estimate, fill = Method)) +
      ggplot2::geom_density(alpha = 0.5) +
      ggplot2::ylab("Density") +
      ggplot2::xlab(ifelse(is.null(y_series_name) | is.null(x_series_name),"Effect Size of the Link Between \n X series and Y series in Daily Life",
                  paste0("Effect Size of the Link Between \n",x_series_name, " and ",y_series_name, " in Daily Life"))) +
      ggplot2::scale_x_continuous(
        breaks = seq(lims[1], lims[2], by = 0.20),  # X axis breaks.
        limits = lims  # Set x-axis limits
      ) +
      ggplot2::theme_bw()

    return(density_plot)

  }

  }

  ### PLOTS FOR LAGGED RELATIONSHIPS.

  if(iarimax_object$type == 'lagged_predictor'){

    #Create Caterpillar Plot.
    if (plot_type == 'caterpillar') {
      # Reorder data
      df_plt <- df_plt %>%
        dplyr::mutate(!!id_var_sym := forcats::fct_reorder(as.factor(!!id_var_sym), xreg_lag))

      # Create color coding
      df_plt <- df_plt %>%
        dplyr::mutate(line_color = dplyr::case_when(
          xreg_lag - crit_val * stderr_xreg_lag > 0 ~ "green",  # Interval entirely positive
          xreg_lag + crit_val * stderr_xreg_lag < 0 ~ "red",    # Interval entirely negative
          TRUE ~ "black"                                # Interval crosses zero or includes it
        ))

      # Create the plot
      caterpillar_plot <- ggplot2::ggplot(df_plt, ggplot2::aes(x = !!id_var_sym, y = xreg_lag)) +
        ggplot2::geom_point() +  # Plot the mean values as points
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::geom_linerange(ggplot2::aes(ymin = xreg_lag - crit_val * stderr_xreg_lag, ymax = xreg_lag + crit_val * stderr_xreg_lag, color = line_color)) +  # Add error bars for +/- 2SD with color based on the new variable
        ggplot2::scale_color_identity() + # Use the colors specified in the iarimax_object
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +  # Add a horizontal line at y = 0 for reference
        ggplot2::coord_flip() +  # Flip the coordinates to make it a horizontal plot
        ggplot2::labs(x = "Participant ID",
                      y = "i-ARIMAX Effect Sizes and 95% Confidence Intervals of \nWithin-Person Associations",
                      title = ifelse(is.null(y_series_name) | is.null(x_series_name), 'X series (t-1) Linked to Y series', paste0(x_series_name, " Linked to  \n", y_series_name)),
                      caption = "Note: Blue vertical lines represent the RE-MA pooled (nomothetic) effect in the middle\n and the lower and upper bounds of the 95% CI of the pooled effect.\nGreen horizontal lines indicate 95% CI of positive associations and red indicate negative.") +
        ggplot2::geom_hline(yintercept = iarimax_object$meta_analysis$beta, color = "blue", linewidth = 1) +  # Add population-level mean intercept line
        ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = iarimax_object$meta_analysis$ci.lb, ymax = iarimax_object$meta_analysis$ci.ub, alpha = 0.2, color = "blue") + # Add 95% credible interval
        ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0)) +
        ggplot2::theme(axis.text.y = ggplot2::element_blank())

      return(caterpillar_plot)

    }

    #Plot empirical cumulative distribution function.
    if (plot_type == "ecdf") {

      if('rand_df' %in% names(iarimax_object)) {
        #Merge random effects to I-arimax dataframe.
        df_plt <- merge(df_plt, subset(iarimax_object$rand_df, select = c(id_var,'random_slope_lag')), by = id_var, all.x=TRUE)
      }
      else {
        stop('Random effects are requiered for ecdf plot. Please run IARIMAXoid_Pro_Lag with hlm_compare = TRUE.')
      }


      #Pivot dataset to long format.
      df_plt_long <- df_plt %>%
        dplyr::select(c(!!id_var_sym,"xreg_lag","random_slope_lag","raw_correlation")) %>%
        dplyr::rename(IarimaxLag = xreg_lag, MLMLag = random_slope_lag, RawCorrLag = raw_correlation) %>%
        tidyr::pivot_longer(
          cols = c("IarimaxLag","MLMLag","RawCorrLag"),
          names_to = "Method",
          values_to = "Estimate"
        )

      #Create plot.
      ecdf_plot <- ggplot2::ggplot(df_plt_long, ggplot2::aes(Estimate, colour = Method)) +
        ggplot2::stat_ecdf(geom = "step", pad = FALSE, linewidth = 1) +
        ggplot2::ylab("Empirical Cumulative Distribution Function") +
        ggplot2::xlab(ifelse(is.null(y_series_name) | is.null(x_series_name),"Effect Size of the Link Between \n X series (t-1) and Y series in Daily Life",
                             paste0("Effect Size of the Link Between \n",x_series_name, " and ",y_series_name, " in Daily Life"))) +
        ggplot2::scale_y_continuous(breaks = seq(0.00, 1.00, by = 0.10)) +
        ggplot2::scale_x_continuous(breaks = seq(-1.0, 1.0, by = 0.20)) +
        ggplot2::theme_bw()

      return(ecdf_plot)

    }

    if (plot_type == "density") {

      if('rand_df' %in% names(iarimax_object)) {
        #Merge random effects to I-arimax dataframe.
        df_plt <- merge(df_plt, subset(iarimax_object$rand_df, select = c(id_var,'random_slope_lag')), by = id_var, all.x=TRUE)
      }
      else {
        stop('Random effects are requiered for density plot. Please run IARIMAXoid_Pro_Lag with hlm_compare = TRUE.')
      }

      #Pivot dataset to long format.
      df_plt_long <- df_plt %>%
        dplyr::select(c(!!id_var_sym,"xreg_lag","random_slope_lag","raw_correlation")) %>%
        dplyr::rename(IarimaxLag = xreg_lag, MLMLag = random_slope_lag, RawCorrLag = raw_correlation) %>%
        tidyr::pivot_longer(
          cols = c("IarimaxLag","MLMLag","RawCorrLag"),
          names_to = "Method",
          values_to = "Estimate"
        )

      #Plot density plot.
      density_plot <- df_plt_long %>%
        ggplot2::ggplot(ggplot2::aes(x = Estimate, fill = Method)) +
        ggplot2::geom_density(alpha = 0.5) +
        ggplot2::ylab("Density") +
        ggplot2::xlab(ifelse(is.null(y_series_name) | is.null(x_series_name),"Effect Size of the Link Between \n X series (t-1) and Y series in Daily Life",
                             paste0("Effect Size of the Link Between \n",x_series_name, " and ",y_series_name, " in Daily Life"))) +
        ggplot2::scale_x_continuous(
          breaks = seq(lims[1], lims[2], by = 0.20),  # X axis breaks.
          limits = lims  # Set x-axis limits
        ) +
        ggplot2::theme_bw()

      return(density_plot)

    }

  }


}

utils::globalVariables(c('line_color','n_valid','n_params','xreg','stderr_xreg',
                         'random_slope','raw_correlation','Estimate','Method',
                         'xreg_lag','stderr_xreg_lag','random_slope_lag')) #Declare global variables.
