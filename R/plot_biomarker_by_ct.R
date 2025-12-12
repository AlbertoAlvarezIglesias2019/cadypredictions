#' Plot Biomarker Trajectory Stratified by Cardiotoxicity (CT) Status
#'
#' Generates a scatter plot and overlays a **LOESS smooth curve** for each
#' **Cardiotoxicity (CT) status** group to visualize the trajectory of a
#' specified biomarker over time (\code{days_from_baseline}).
#' This function is designed to compare biomarker changes between
#' participants who developed CT ("Yes") and those who did not ("No"), based on the
#' binary grouping variable specified by \code{ct_col}.
#'
#' @param data A data frame containing the data. It must include the columns
#'   \code{days_from_baseline} (numeric, representing time),
#'   the column specified by \code{biomarker_name}, and the grouping column
#'   specified by \code{ct_col}.
#' @param biomarker_name A character string specifying the name of the column
#'   containing the numeric biomarker values to be plotted on the y-axis (e.g.,
#'   \code{"Troponin"}, \code{"NTproBNP_change"}).
#' @param ct_col A character string specifying the name of the column containing
#'   the binary grouping variable for Cardiotoxicity (CT) status.
#'   This column is expected to have levels **"Yes"** and **"No"**.
#'   Defaults to \code{"CT_mp_50_YN"}.
#' @param x_lab A character string for the x-axis label. Defaults to
#'   \code{"Days from First Sample"}.
#' @param y_lab A character string for the y-axis label. Defaults to
#'   \code{"Change in Biomarker"}.
#' @param title A character string for the main title of the plot. If \code{NULL} (default),
#'   a title is automatically generated using \code{biomarker_name}.
#'
#' @return A \code{ggplot} object displaying the biomarker's trajectory over time,
#'   with separate LOESS smoothers and coloring for the two CT status groups.
#'
#' @details
#' \itemize{
#'   \item **LOESS Smoother**: LOESS (Locally Estimated Scatterplot Smoothing)
#'         lines are generated separately for the "Yes" (developed CT) and "No"
#'         (did not develop CT) groups.
#'   \item **Aesthetics**: The plot uses manually defined colors for consistency:
#'         "No" (no CT) is mapped to a **cyan/blue** color (`#00BFC4`), and "Yes"
#'         (developed CT) is mapped to a **red** color (`#F8766D`). 
#'   \item **Non-Standard Evaluation (NSE)**: Both \code{biomarker_name} and
#'         \code{ct_col} arguments are used via \code{.data[[column_name]]} to select
#'         the columns for the y-axis, color, and fill, allowing the column names
#'         to be passed as strings.
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth labs theme_minimal theme
#'   element_text element_rect element_blank scale_colour_manual scale_fill_manual
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' # Example 1: Using the default CT column name
#' # Assuming 'cardio_data' has columns: 'days_from_baseline', 'NTproBNP', 'CT_mp_50_YN'
#'
#' plot_biomarker_by_ct(
#'   data = cardio_data,
#'   biomarker_name = "NTproBNP",
#'   y_lab = "NT-proBNP Concentration (ng/mL)"
#' )
#'
#' # Example 2: Specifying a custom CT column name
#' # Assuming 'cardio_data' also has a column 'CT_mp_30_YN'
#'
#' plot_biomarker_by_ct(
#'   data = cardio_data,
#'   biomarker_name = "Troponin",
#'   ct_col = "CT_mp_30_YN",
#'   y_lab = "Troponin Concentration (ng/mL)",
#'   title = "Troponin Trajectory by CT (30% MP) Status"
#' )
#' }
#'
plot_biomarker_by_ct <- function(data,
                                 biomarker_name,
                                 ct_col = "CT_mp_50_YN",
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
      colour = .data[[ct_col]],
      fill   = .data[[ct_col]]
    )
  ) +
    geom_point(alpha = 0.1, size = 1.8) +
    geom_smooth(method = "loess", se = TRUE, alpha = 0.33) +
    scale_colour_manual(
      values = c("No" = "#00BFC4",
                 "Yes" = "#F8766D")
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