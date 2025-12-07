#' Plot the Overall Trajectory of Biomarker Change Over Time
#'
#' Generates a scatter plot and overlays a **LOESS smooth curve** to visualize
#' the change (\code{delta}) in a specified biomarker over time (\code{days_from_baseline})
#' across all treatment groups combined. This function focuses on the **general trend**
#' of the biomarker's trajectory rather than group comparisons or individual trajectories.
#'
#' @param data A data frame in long format. It must contain the columns \code{biomarker},
#'   \code{days_from_baseline} (numeric, representing time), and \code{delta} (numeric,
#'   representing the change from baseline or a reference point).
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
#'   biomarker over time with an overlaid LOESS smoother.
#'
#' @details
#' \itemize{
#'   \item **Data Filtering**: The function first filters the input \code{data}
#'         to include only rows matching the specified \code{biomarker_name}.
#'   \item **LOESS Smoother**: A **red** LOESS (Locally Estimated Scatterplot Smoothing)
#'         line is used to depict the non-linear trend, including the standard
#'         error (SE) band. 
#'   \item **Scatter Points**: Individual data points are displayed in **steelblue**
#'         with high transparency (\code{alpha = 0.1}) to handle potentially large
#'         datasets without visual clutter.
#'   \item **Aesthetics**: Uses \code{theme_minimal} and centers the legend in the
#'         top-right corner.
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth labs theme_minimal theme
#'   element_text element_rect element_blank filter
#'
#' @examples
#' \dontrun{
#' # Assuming 'long_biomarker_data' has the required columns:
#' # 'biomarker', 'days_from_baseline', 'delta'
#'
#' plot_biomarker_delta_overall(
#'   data = long_biomarker_data,
#'   biomarker_name = "IL-6",
#'   y_lab = "IL-6 Concentration Change (pg/mL)",
#'   title = "Overall Trajectory of IL-6 Change"
#' )
#' }
#'
plot_biomarker_delta_overall <- function(data,
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
             y = delta)) +
    geom_point(color = "steelblue", alpha = 0.1, size = 2) +
    geom_smooth(method = "loess", se = TRUE, color = "#CD0000", size = 1) +
    labs(
      title = title,
      x = x_lab,
      y = y_lab
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