#' Plot General Distribution of a Numeric Variable by Group
#'
#' This function generates a **frequency histogram** for a specified numeric
#' column (\code{value}), faceted by the required \code{Group} column. It
#' optionally filters the data based on a secondary \code{Variable} column.
#' Vertical dashed lines are added to indicate the **mean** and **median**
#' values of the \code{value} within each group.
#'
#' @param data A data frame containing the data to be plotted. It must include
#'   a column named \code{Group} for faceting and the column specified by \code{value}.
#'   If \code{variable} is used, it must also contain a column named \code{Variable}.
#' @param variable A character string specifying a value in the \code{Variable}
#'   column to filter the data by (e.g., "ALT", "Creatinine"). If \code{NULL} (default),
#'   no filtering is applied based on a \code{Variable} column.
#' @param value A **bare column name** (non-standard evaluation) representing the
#'   numeric variable to be plotted on the x-axis (e.g., \code{Value}, \code{Dose}).
#'   This column is expected to be convertible to numeric.
#' @param title A character string for the main title of the generated plot.
#' @param binwidth A numeric value specifying the width of the bins to use for
#'   the histogram component of the plot.
#'
#' @return A \code{ggplot} object displaying the frequency distribution of the
#'   specified numeric variable, faceted by the \code{Group} column.
#'
#' @details
#' The function relies on the **tidyverse** (specifically \code{dplyr} and
#' \code{ggplot2}) for data manipulation and plotting.
#'
#' \itemize{
#'   \item The y-axis displays **Frequency** (\code{after_stat(count)}) rather than density.
#'   \item The \code{value} argument uses **non-standard evaluation (NSE)**; it should
#'         be supplied as a bare column name (e.g., \code{Value}).
#'   \item Summary statistics (mean and median) are calculated and plotted for each \code{Group}.
#'   \item Mean is plotted as a \strong{blue dashed line}, and median as a \strong{red dashed line}.
#' }
#'
#' @importFrom dplyr filter group_by summarize %>%
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline facet_wrap
#'   theme_minimal labs scale_color_manual after_stat
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' # Example 1: Plotting a generic 'Value' column across all data
#' # assuming 'my_data' has 'Group' and 'Value' columns
#' plot_general_distribution(
#'   data = my_data,
#'   value = Value,
#'   title = "Overall Value Distribution by Group",
#'   binwidth = 10
#' )
#'
#' # Example 2: Filtering data for a specific 'Variable' ("Sodium")
#' # assuming 'long_data' has 'Group', 'Variable', and 'Measurement' columns
#' plot_general_distribution(
#'   data = long_data,
#'   variable = "Sodium",
#'   value = Measurement,
#'   title = "Sodium Distribution by Group",
#'   binwidth = 2
#' )
#' }
#'

plot_general_distribution <- function(data, variable = NULL, value, title, binwidth) {
  
  # Ensure 'value' is treated as a column name (string)
  value <- deparse(substitute(value))
  
  # If variable is provided, filter the dataset
  if (!is.null(variable)) {
    data <- data %>% filter(Variable == variable)
  }
  
  # Compute summary statistics
  disp_summary <- data %>%
    group_by(Group) %>%
    summarize(mean_val = mean(as.numeric(.data[[value]]), na.rm = TRUE),
              median_val = median(as.numeric(.data[[value]]), na.rm = TRUE), .groups = "drop")
  
  # Generate the histogram and density plot
  plot <- data %>%
    ggplot(aes(x = as.numeric(.data[[value]]))) +
    geom_histogram(aes(y = after_stat(count)), binwidth = binwidth, fill = "skyblue", color = "black", alpha = 0.5) +
    # geom_density(aes(y = after_stat(count)), color = "black", size = 0.5) +
    geom_vline(data = disp_summary, aes(xintercept = mean_val, color = "Mean"), linetype = "dashed", size = 1) +
    geom_vline(data = disp_summary, aes(xintercept = median_val, color = "Median"), linetype = "dashed", size = 1) +
    facet_wrap(~ Group) +
    theme_minimal() +
    labs(title = title, x = ifelse(is.null(variable), value, variable), y = "Frequency") +
    scale_color_manual(values = c("Mean" = "blue", "Median" = "red"))
  
  return(plot)
}