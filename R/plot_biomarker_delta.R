#' Plot the Trajectory of Biomarker Change, Stratified by Treatment Regimen
#'
#' Generates a scatter plot and overlays a **LOESS smooth curve** for each treatment
#' group to visualize the change (\code{delta}) in a specified biomarker over time
#' (\code{days_from_baseline}). This function allows for direct comparison of the
#' biomarker trajectory between different treatment regimens.
#'
#' @param data A data frame in long format. It must contain the columns \code{biomarker},
#'   \code{days_from_baseline} (numeric, representing time), \code{delta} (numeric,
#'   representing the change from baseline or a reference point), and \code{treatment_reg}
#'   (the categorical grouping variable).
#' @param biomarker_name A character string specifying the specific value in the
#'   \code{biomarker} column to plot (e.g., "CRP", "Ferritin").
#' @param x_lab A character string for the x-axis label. Defaults to
#'   \code{"Days from First Sample"}.
#' @param y_lab A character string for the y-axis label. Defaults to
#'   \code{"Change in Biomarker"}.
#' @param title A character string for the main title of the plot. If \code{NULL} (default),
#'   a title is automatically generated using \code{biomarker_name}.
#'
#' @return A \code{ggplot} object displaying the scatter plot of the change in the
#'   biomarker over time, with separate LOESS smoothers and coloring for each
#'   \code{treatment_reg} group.
#'
#' @details
#' \itemize{
#'   \item **Data Filtering**: The function first filters the input \code{data}
#'         to include only rows matching the specified \code{biomarker_name}.
#'   \item **Grouped LOESS Smoother**: Both the scatter points and the LOESS
#'         smooth curves are colored and grouped by the \code{treatment_reg} variable.
#'         The LOESS smoothers are calculated separately for each group, allowing
#'         for visual assessment of the differential trajectory over time. 
#'   \item **Aesthetics**: Uses \code{theme_minimal} and places the legend inside
#'         the top-right corner. The scatter points use low transparency (\code{alpha = 0.1}).
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth labs theme_minimal theme
#'   element_text element_rect element_blank filter
#'
#' @examples
#' \dontrun{
#' # Assuming 'long_biomarker_data' has the required columns:
#' # 'biomarker', 'days_from_baseline', 'delta', 'treatment_reg'
#'
#' plot_biomarker_delta(
#'   data = long_biomarker_data,
#'   biomarker_name = "Interleukin-6",
#'   y_lab = "IL-6 Concentration Change (pg/mL)",
#'   title = "IL-6 Change Over Time by Treatment Regimen"
#' )
#' }
#'
plot_biomarker_delta <- function(data,
                                 biomarker_name,
                                 x_lab = "Days from First Sample",
                                 y_lab = "Change in Biomarker",
                                 title = NULL) {
  
  # Default title if not provided
  if (is.null(title)) {
    title <- paste("Change in", biomarker_name, "over Time by Treatment Group")
  }
  
  ggplot(data %>% filter(biomarker == biomarker_name),
         aes(x = days_from_baseline,
             y = delta,
             color = treatment_reg)) +
    geom_point(alpha = 0.1, size = 2) +
    geom_smooth(method = "loess", se = TRUE, aes(group = treatment_reg)) +
    labs(
      title = title,
      x = x_lab,
      y = y_lab,
      color = ""
    ) +
    theme_minimal(base_size = 14) +
    theme(
      #plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      legend.position = c(0.90, 0.90),    # legend inside top-right corner
      # legend.background = element_rect(fill = alpha("white", 0.6), color = NA),
      legend.key = element_rect(fill = NA),
      legend.text = element_text(size = 12)
    )
}