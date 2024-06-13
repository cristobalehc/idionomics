#' Create idionomic plots.
#'
#' @export
#' @importFrom utils globalVariables
#'
#' @param dataframe Your dataframe.
#' @param id_var Your id variable.
#' @param y_series_name Optional: A string containing the substantive name of your y series to be used in the plot title.
#' @param x_series_name Optional: A string containing the substantive name of your x series to be used in the plot title.
#' @param alpha_crit_t Critical value for model significance. Defaluts to 0.05.
#' @param plot_type Plot type. Currently supported: "caterpillar".
#'
#' @returns A list containing a dataframe with the ARIMA parameteres, plus the xreg parameter (the beta value for your x_series) together with their std.errors. If metaanalysis = TRUE, will also output a random effects meta analysis. If hlm_compare = TRUE, will also output a model comparison with HLM.


#######################################
############ PLOTOMATRON #######
#####################################



Plotomatron <- function(dataframe,id_var ,y_series_name = NULL, x_series_name = NULL, alpha_crit_t = 0.05, plot_type = 'caterpillar') {

  # Check if df is an iarimax object
  if (is.list(dataframe) && all(c('results_df', 'meta_analysis') %in% names(dataframe))) {
    df_plt <- dataframe$results_df
  } else {
    stop('Plotomatron requieres an iarimaxoid object that contains "results_df" and "meta_analysis". Run IARIMAXoid_Pro with metaanalysis = TRUE')
  }


    #Create symbolic variables.
    xreg_sym <- rlang::sym('xreg')
    xregstd_sym <- rlang::sym('stderr_xreg')
    n_valid_sym <- rlang::sym('n_valid')
    n_params_sym <- rlang::sym('n_params')
    id_var_sym <- rlang::sym(id_var)


  #Calcualte degrees of freedom and critical value.
  df_plt <- df_plt %>%
    dplyr::mutate(df_mod = !!n_valid_sym - !!n_params_sym)
  #Critical value.
  crit_val <- stats::qt(1 - (alpha_crit_t / 2), df_plt$df_mod)

  if (plot_type == 'caterpillar') {
    # Reorder data
    df_plt <- df_plt %>%
      dplyr::mutate(!!id_var_sym := forcats::fct_reorder(as.factor(!!id_var_sym), !!xreg_sym))

    # Create color coding
    df_plt <- df_plt %>%
      dplyr::mutate(line_color = dplyr::case_when(
        !!xreg_sym - crit_val * !!xregstd_sym > 0 ~ "green",  # Interval entirely positive
        !!xreg_sym + crit_val * !!xregstd_sym < 0 ~ "red",    # Interval entirely negative
        TRUE ~ "black"                                # Interval crosses zero or includes it
      ))

    # Create the plot
    ggplot2::ggplot(df_plt, ggplot2::aes(x = !!id_var_sym, y = !!xreg_sym)) +
      ggplot2::geom_point() +  # Plot the mean intercepts as points
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::geom_linerange(ggplot2::aes(ymin = !!xreg_sym - crit_val * !!xregstd_sym, ymax = !!xreg_sym + crit_val * !!xregstd_sym, color = line_color)) +  # Add error bars for +/- 2SD with color based on the new variable
      ggplot2::scale_color_identity() + # Use the colors specified in the dataframe
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +  # Add a horizontal line at y = 0 for reference
      ggplot2::coord_flip() +  # Flip the coordinates to make it a horizontal plot
      ggplot2::labs(x = "Participant ID",
                    y = "i-ARIMAX Effect Sizes and 95% Confidence Intervals of Within-Person Associations",
                    title = ifelse(is.null(y_series_name) | is.null(x_series_name), 'X series Linked to Y series', paste0(x_series_name, " Linked to  \n", y_series_name)),
                    caption = "Note: Blue vertical lines represent the RE-MA pooled (nomothetic) effect in the middle and the lower and upper bounds of the 95% CI of the pooled effect.\nGreen horizontal lines indicate 95% CI of positive associations and red indicate negative ones.") +
      ggplot2::geom_hline(yintercept = dataframe$meta_analysis$beta, color = "blue", linewidth = 1) +  # Add population-level mean intercept line
      ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = dataframe$meta_analysis$ci.lb, ymax = dataframe$meta_analysis$ci.ub, alpha = 0.2, color = "blue") + # Add 95% credible interval
      ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0)) +
      ggplot2::theme(axis.text.y = ggplot2::element_blank())

  }
  else {
    stop('plot_type = "caterpillar" is currently the only option available')
  }
}

utils::globalVariables(c('line_color')) #Declare global variables.
