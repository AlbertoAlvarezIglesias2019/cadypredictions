#' @title Generate a LOESS smoothed plot for longitudinal data
#'
#' @description This function creates a \code{ggplot2} line plot, typically for longitudinal data,
#'   with lines grouped by \code{group_var} and colored by \code{color_group}. It overlays a
#'   LOESS (Locally Estimated Scatterplot Smoothing) curve for each color group to visualize
#'   trends.
#'
#' @param data A data frame containing the data to be plotted.
#' @param x_var A character string specifying the name of the column to be used for the x-axis (e.g., time).
#' @param y_var A character string specifying the name of the column to be used for the y-axis (e.g., a measurement).
#' @param group_var A character string specifying the column that defines individual lines (e.g., SubjectID).
#' @param color_group A character string specifying the column used to color the LOESS smooths and lines.
#'   This variable is treated as a factor.
#' @param plot_title A character string for the main title of the plot.
#' @param x_label A character string for the x-axis label.
#' @param y_label A character string for the y-axis label.
#' @param x_limits An optional numeric vector of length two (e.g., \code{c(min, max)}) to set the
#'   limits of the x-axis using \code{xlim()}. Defaults to \code{NULL}.
#' @param y_limits An optional numeric vector of length two (e.g., \code{c(min, max)}) to set the
#'   limits of the y-axis using \code{ylim()}. Defaults to \code{NULL}.
#'
#' @return The function prints the generated \code{ggplot} object but does not explicitly return it
#'   for further manipulation.
#'
#' @details
#' \itemize{
#'   \item Individual subject lines (\code{group_var}) are plotted in a faint gray color (\code{alpha = 0.4}).
#'   \item The primary visual focus is the LOESS smoothed curve (\code{geom_smooth(method = "loess")}),
#'     which includes standard error (SE) ribbons.
#'   \item The plot uses a \code{theme_minimal()} and places the legend at the bottom.
#' }
#'
#' @importFrom ggplot2 ggplot aes sym geom_line geom_smooth theme_minimal labs theme element_blank guides guide_legend xlim ylim
#' @importFrom rlang !!
#'
#' @examples
#' \dontrun{
#' # Assuming 'my_data' is a dataframe with columns: SubjectID, Visit, Measurement, Treatment
#' loess_plot(
#'   data = my_data,
#'   x_var = "Visit",
#'   y_var = "Measurement",
#'   group_var = "SubjectID",
#'   color_group = "Treatment",
#'   plot_title = "Longitudinal Trends by Treatment Group",
#'   x_label = "Study Day",
#'   y_label = "Biomarker Level",
#'   x_limits = c(0, 100)
#' )
#' }
#' @export


loess_plot <- function(data, x_var, y_var, group_var, color_group, plot_title, x_label, y_label,
                       x_limits = NULL, y_limits = NULL) {
  x_sym <- sym(x_var)
  y_sym <- sym(y_var)
  group_sym <- sym(group_var)
  color_sym <- sym(color_group)

  plot <- ggplot(data, aes(x = !!x_sym, y = !!y_sym, group = !!group_sym, color = as.factor(!!color_sym))) +
    geom_line(color = "lightgray", alpha = 0.4) +
    geom_smooth(method = "loess", aes(group = as.factor(!!color_sym), color = as.factor(!!color_sym)), se = TRUE) +
    theme_minimal() +
    labs(x = x_label, y = y_label, title = plot_title) +
    theme(
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "bottom"  # Move legend to the bottom
    ) +
    guides(color = guide_legend(title = NULL))

  # Apply x and y limits only if provided
  if (!is.null(x_limits)) {
    plot <- plot + xlim(x_limits)
  }
  if (!is.null(y_limits)) {
    plot <- plot + ylim(y_limits)
  }

  print(plot)
}
