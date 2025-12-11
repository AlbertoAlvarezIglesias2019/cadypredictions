#' Create a Set of Scatter Plots for a Predictor Against Multiple Biomarkers
#'
#' Generates a named list of five individual scatter plots, each comparing a
#' single \code{predictor} variable against a fixed set of five baseline
#' biomarkers. This function relies on the helper functions \code{scatter_plot()}
#' and \code{format_corr()} (via \code{corr_list}) to create and annotate each plot.
#'
#' @param data A data frame containing the \code{predictor} variable and all five
#'   biomarker columns.
#' @param predictor A character string specifying the name of the column to be
#'   used as the predictor variable for the x-axis in all plots (e.g., "Age", "BMI").
#' @param corr_list A named list of five character strings, where the names match
#'   the biomarker columns and the values are the formatted correlation results
#'   for that biomarker vs. the predictor (e.g., the output of \code{get_corr_strings()}).
#'
#' @return A **named list** of five \code{ggplot} objects. The names of the list
#'   elements correspond to the biomarker names: \code{"BNP_bl"}, \code{"NT_pro_BNP_bl"},
#'   \code{"CRP_bl"}, \code{"hsTnI_STAT_bl"}, and \code{"Galectin_3_bl"}.
#'
#' @details
#' \itemize{
#'   \item **Biomarkers**: The function is designed to plot against a fixed set of
#'     five biomarker columns, which must match the names in \code{corr_list}:
#'     \code{"BNP_bl"}, \code{"NT_pro_BNP_bl"}, \code{"CRP_bl"}, \code{"hsTnI_STAT_bl"},
#'     and \code{"Galectin_3_bl"}.
#'   \item **Plot Generation**: It uses \code{mapply} to iterate over the biomarker
#'     names, custom y-axis labels, and specific y-axis limits, calling the
#'     \code{\link{scatter_plot}} helper function for each combination.
#'   \item **Y-Limits**: Custom y-axis limits are applied to the "hs-Troponin I"
#'     plot (\code{c(0, 15)}) while other plots use the default limits (\code{NULL}).
#'   \item **Annotation**: The formatted correlation result string from \code{corr_list}
#'     is used as the **subtitle** for each plot, providing the statistical summary.
#' }
#'
#' @seealso
#' \code{\link{scatter_plot}} for plot details, \code{\link{format_corr}} for string formatting, and \code{\link{get_corr_strings}} for creating \code{corr_list}.
#'
#' @importFrom purrr map set_names
#' @importFrom rlang all_of
#'
#' @examples
#' \dontrun{
#' # Assuming 'clinical_data' and 'age_corr_strings' (from get_corr_strings example) are available
#'
#' # Create a list of 5 scatter plots comparing 'Age' to the biomarkers
#' age_plot_grid <- build_plots(
#'   data = clinical_data,
#'   predictor = "Age",
#'   corr_list = age_corr_strings
#' )
#'
#' # The list can then be arranged into a grid (e.g., using patchwork::wrap_plots)
#' # patchwork::wrap_plots(age_plot_grid, ncol = 3)
#' }
#'
build_plots <- function(data, predictor, corr_list) {
  
  biomks  <- names(corr_list)
  y_labels <- c("Baseline BNP", "Baseline NT-proBNP", "Baseline CRP", "Baseline hs-Troponin I", "Baseline Galectin-3")
  
  ylim_opt <- list(NULL, NULL, NULL, c(0, 15), NULL)
  
  plots <- mapply(
    function(b, lab, yl)
      scatter_plot(data, predictor, b, corr_list[[b]], lab, yl),
    biomks, y_labels, ylim_opt,
    SIMPLIFY = FALSE
  )
  
  # Return as named list
  names(plots) <- biomks
  plots
}