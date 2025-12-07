#' Plot Biomarker Trajectory Stratified by Cardiotoxicity (CT) Status
#'
#' Generates a scatter plot and overlays a **LOESS smooth curve** for each
#' **Cardiotoxicity (CT) status** group to visualize the trajectory of a
#' specified biomarker over time (\code{days_from_baseline}).
#' This function is specifically designed to compare biomarker changes between
#' participants who developed CT (Yes) and those who did not (No), based on the
#' binary variable \code{CT_mp_50_YN}.
#'
#' @param data A data frame containing the data. It must include the columns
#'   \code{days_from_baseline} (numeric, representing time),
#'   \code{CT_mp_50_YN} (the binary factor/character grouping variable, expected
#'   to have levels "Yes" and "No"), and the column specified by \code{biomarker_name}.
#' @param biomarker_name A character string specifying the name of the column
#'   containing the numeric biomarker values to be plotted on the y-axis (e.g.,
#'   \code{"Troponin"}, \code{"NTproBNP_change"}).
#' @param x_lab A character string for the x-axis label. Defaults to
#'   \code{"Days from First Sample"}.
#' @param y_lab A character string for the y-axis label. Defaults to
#'   \code{"Change in Biomarker"}.
#' @param title A character string for the main title of the plot. If \code{NULL} (default),
#'   a title is automatically generated using \code{biomarker_name}.
#'
#' @return A \code{ggplot} object displaying the biomarker's trajectory over time,
#'   with separate LOESS smoothers and coloring for the two CT status groups.
#'
#' @details
#' \itemize{
#'   \item **LOESS Smoother**: LOESS (Locally Estimated Scatterplot Smoothing)
#'         lines are generated separately for the "Yes" (developed CT) and "No"
#'         (did not develop CT) groups.
#'   \item **Aesthetics**: The plot uses manually defined colors for consistency:
#'         "No" (no CT) is mapped to a **cyan/blue** color (`#00BFC4`), and "Yes"
#'         (developed CT) is mapped to a **red** color (`#F8766D`). 
#'   \item **Non-Standard Evaluation (NSE)**: The \code{biomarker_name} argument is
#'         used via \code{.data[[biomarker_name]]} to select the column for the y-axis,
#'         allowing the column name to be passed as a string.
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth labs theme_minimal theme
#'   element_text element_rect element_blank scale_colour_manual scale_fill_manual
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' # Assuming 'cardio_data' has columns: 'days_from_baseline', 'NTproBNP', 'CT_mp_50_YN'
#'
#' plot_biomarker_by_ct(
#'   data = cardio_data,
#'   biomarker_name = "NTproBNP",
#'   y_lab = "NT-proBNP Concentration (ng/mL)",
#'   title = "NT-proBNP Trajectory by Cardiotoxicity Status"
#' )
#' }
#'
plot_biomarker_by_ct <- function(data,
                                 biomarker_name,
                                 x_lab = "Days from First Sample",
                                 y_lab = "Change in Biomarker",
                                 title = NULL) {
  
  # Default title if not provided
  if (is.null(title)) {
    title <- paste("Change in", biomarker_name, "over Time by CT status")
  }
  
  ggplot(
    data,
    aes(
      x = days_from_baseline,
      y = .data[[biomarker_name]],
      colour = CT_mp_50_YN,
      fill  = CT_mp_50_YN
    )
  ) +
    geom_point(alpha = 0.2, size = 1.8) +  # keep mapped colours
    geom_smooth(method = "loess", se = TRUE, alpha = 0.5) + # Changed alfa to alpha
    scale_colour_manual(
      values = c("No" = "#00BFC4",  # default ggplot blue
                 "Yes" = "#F8766D") # default ggplot red
    ) +
    scale_fill_manual(
      values = c("No" = "#00BFC4",
                 "Yes" = "#F8766D")
    ) +
    labs(
      title = title,
      x = x_lab,
      y = y_lab,
      colour = "CT status",
      fill = "CT status"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = c(0.90, 0.90),
      legend.key = element_rect(fill = NA),
      legend.text = element_text(size = 12)
    )
}