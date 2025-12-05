#' Plot the Distribution of a Laboratory Variable by Treatment Arm
#'
#' Generates a histogram and density plot for a specified laboratory variable
#' at a particular study visit, faceted by the Randomized Arm (\code{RandArm}).
#' The plot includes vertical dashed lines indicating the mean and median
#' values for the variable within each randomized arm.
#'
#' @param data A tidy data frame containing at least the columns \code{Visit},
#'   \code{Variable}, \code{Value}, and \code{RandArm}.
#' @param lab_variable A character string specifying the value in the \code{Variable}
#'   column to plot (e.g., "ALT", "Creatinine"). This will also be used as the x-axis label.
#' @param visit A character string or factor specifying the value in the \code{Visit}
#'   column to filter the data by (e.g., "Week 4", "Baseline").
#' @param title A character string for the main title of the generated plot.
#' @param binwidth A numeric value specifying the width of the bins to use for
#'   the histogram component of the plot.
#'
#' @return A \code{ggplot} object displaying the distribution of the specified
#'   laboratory variable, faceted by \code{RandArm}.
#'
#' @details
#' This function requires the \strong{tidyverse} package (specifically \code{dplyr}
#' and \code{ggplot2}).
#' The mean is represented by a \strong{blue dashed line} and the median by a
#' \strong{red dashed line}.
#' The y-axis of the histogram is scaled to density (\code{..density..}) to align
#' with the overlaid kernel density estimate (\code{geom_density}).
#'
#' @importFrom dplyr filter group_by summarize %>%
#' @importFrom ggplot2 ggplot aes geom_histogram geom_density geom_vline facet_wrap
#'   theme_minimal labs scale_color_manual
#'
#' @examples
#' \dontrun{
#' # Assuming 'clinical_data' is a suitable data frame
#' # Plotting Creatinine distribution at Week 12 with a binwidth of 5
#' creatinine_plot <- plot_distribution(
#'   data = clinical_data,
#'   lab_variable = "Creatinine",
#'   visit = "Week 12",
#'   title = "Creatinine Distribution at Week 12 by Treatment Arm",
#'   binwidth = 5
#' )
#' print(creatinine_plot)
#' }
#'
plot_distribution <- function(data, lab_variable, visit, title, binwidth) {
  # Filter for the specific visit and lab variable
  disp_summary <- data %>%
    filter(Visit == visit, Variable == lab_variable) %>%
    group_by(RandArm) %>%
    summarize(mean_val = mean(Value, na.rm = TRUE),
              median_val = median(Value, na.rm = TRUE), .groups = "drop")
  
  # Generate the histogram and density plot
  plot <- data %>%
    filter(Visit == visit, Variable == lab_variable) %>%
    ggplot(aes(x = Value)) +
    geom_histogram(aes(y = ..density..), binwidth = binwidth, fill = "skyblue", color = "black", alpha = 0.5) +
    geom_density(color = "black", size = 0.5) +
    geom_vline(data = disp_summary, aes(xintercept = mean_val, color = "Mean"), linetype = "dashed", size = 1) +
    geom_vline(data = disp_summary, aes(xintercept = median_val, color = "Median"), linetype = "dashed", size = 1) +
    facet_wrap(~ RandArm) +
    theme_minimal() +
    labs(title = title, x = lab_variable, y = "Density") +
    scale_color_manual(values = c("Mean" = "blue", "Median" = "red"))
  
  return(plot)
}



plot_distribution <- function(data, lab_variable, visit, title, binwidth) {
  # Filter for the specific visit and lab variable
  disp_summary <- data %>%
    filter(Visit == visit, Variable == lab_variable) %>%
    group_by(RandArm) %>%
    summarize(mean_val = mean(Value, na.rm = TRUE),
              median_val = median(Value, na.rm = TRUE), .groups = "drop")
  
  # Generate the histogram and density plot
  plot <- data %>%
    filter(Visit == visit, Variable == lab_variable) %>%
    ggplot(aes(x = Value)) +
    geom_histogram(aes(y = ..density..), binwidth = binwidth, fill = "skyblue", color = "black", alpha = 0.5) +
    geom_density(color = "black", size = 0.5) +
    geom_vline(data = disp_summary, aes(xintercept = mean_val, color = "Mean"), linetype = "dashed", size = 1) +
    geom_vline(data = disp_summary, aes(xintercept = median_val, color = "Median"), linetype = "dashed", size = 1) +
    facet_wrap(~ RandArm) +
    theme_minimal() +
    labs(title = title, x = lab_variable, y = "Density") +
    scale_color_manual(values = c("Mean" = "blue", "Median" = "red"))
  
  return(plot)
}

