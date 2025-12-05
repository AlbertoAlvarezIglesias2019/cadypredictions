#' Plot Individual Trajectories with an Overall LOESS Smoother
#'
#' This function generates a plot showing individual participant trajectories
#' for two variables (\code{x_var} and \code{y_var}) over time or another grouping
#' factor. Crucially, it overlays a single, overall LOESS (Locally Estimated
#' Scatterplot Smoothing) curve and its standard error band to show the
#' general trend across all individuals.
#'
#' @param data A data frame containing the data to be plotted.
#' @param x_var A character string specifying the name of the column for the
#'   x-axis (e.g., "Time", "Visit_Number").
#' @param y_var A character string specifying the name of the column for the
#'   y-axis (the measured outcome).
#' @param group_var A character string specifying the name of the column that
#'   identifies unique individuals or groups (e.g., "SubjectID"). This is
#'   used for the individual \code{geom_line} trajectories.
#' @param plot_title A character string for the main title of the plot.
#' @param x_label A character string for the x-axis label.
#' @param y_label A character string for the y-axis label.
#' @param x_limits A numeric vector of length 2 (e.g., \code{c(min, max)}) to
#'   set the limits of the x-axis. Defaults to \code{NULL} (uses \code{ggplot2} defaults).
#' @param y_limits A numeric vector of length 2 (e.g., \code{c(min, max)}) to
#'   set the limits of the y-axis. Defaults to \code{NULL} (uses \code{ggplot2} defaults).
#' @param smooth_color A color string (e.g., \code{"#0072B2"}, \code{"blue"})
#'   for the LOESS smooth line and its standard error band.
#' @param line_alpha A numeric value between 0 and 1 controlling the transparency
#'   of the individual trajectory lines. A low value (e.g., 0.15) is typical for
#'   plots with many individuals.
#' @param smooth_size A numeric value controlling the line thickness of the LOESS
#'   smooth curve.
#' @param span_value A numeric value (between 0 and 1) passed to the \code{span}
#'   argument of \code{geom_smooth(method = "loess")}. This controls the degree
#'   of smoothing; smaller values result in a wigglier, less smooth curve.
#'
#' @return The function prints the generated \code{ggplot} object to the console
#'   and implicitly returns it.
#'
#' @details
#' This function uses **non-standard evaluation (NSE)** with \code{rlang::sym()}
#' to allow variable names to be passed as character strings.
#'
#' The individual participant lines (\code{geom_line}) are drawn with
#' \code{color = "gray60"} and controlled transparency (\code{line_alpha}).
#'
#' The overall LOESS smoother is applied to the entire dataset (\code{aes(group = 1)})
#' and is not faceted or grouped by an arm/treatment variable.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_smooth theme_minimal labs
#'   theme element_blank element_text margin xlim ylim
#' @importFrom rlang sym
#' @importFrom scales alpha
#'
#' @examples
#' \dontrun{
#' # Assuming 'long_data' has columns: 'SubjectID', 'TimePoint', 'Measurement'
#' loess_plot_overall(
#'   data = long_data,
#'   x_var = "TimePoint",
#'   y_var = "Measurement",
#'   group_var = "SubjectID",
#'   plot_title = "Overall LOESS Trend of Measurement Over Time",
#'   x_label = "Study Time (Days)",
#'   y_label = "Observed Measurement (Units)",
#'   smooth_color = "darkgreen",
#'   span_value = 0.7
#' )
#' }
#'
loess_plot_overall <- function(data, x_var, y_var, group_var,
                               plot_title, x_label, y_label,
                               x_limits = NULL, y_limits = NULL,
                               smooth_color = "#0072B2",
                               line_alpha = 0.15,
                               smooth_size = 1.3,
                               span_value = 0.5) {
  
  x_sym <- sym(x_var)
  y_sym <- sym(y_var)
  group_sym <- sym(group_var)
  
  plot <- ggplot(data, aes(x = !!x_sym, y = !!y_sym, group = !!group_sym)) +
    # Participant trajectories (grouped)
    geom_line(color = "gray60", alpha = line_alpha) +
    
    # Global smoother across all data (not grouped)
    geom_smooth(
      aes(group = 1),
      method = "loess", se = TRUE,
      color = smooth_color,
      fill = scales::alpha(smooth_color, 0.25),
      linewidth = smooth_size,
      span = span_value
    ) +
    
    theme_minimal(base_size = 13) +
    labs(
      title = plot_title,
      x = x_label,
      y = y_label
    ) +
    theme(
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10))
    )
  
  if (!is.null(x_limits)) plot <- plot + xlim(x_limits)
  if (!is.null(y_limits)) plot <- plot + ylim(y_limits)
  
  print(plot)
}