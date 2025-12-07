#' Create a Box Plot with Optional Summary Table
#'
#' Generates a box plot for a specified numeric variable (\code{y_variable}),
#' grouped by a categorical variable (\code{group}). It optionally filters the data
#' by a single \code{visit} and can automatically combine the plot with a table
#' of summary statistics (n, mean, sd, median, IQR) using the \code{ggpubr} package.
#'
#' @param data A data frame containing the data to be plotted.
#' @param group A **bare column name** (non-standard evaluation) representing the
#'   categorical variable to be plotted on the x-axis and used for grouping
#'   (e.g., \code{RandArm}, \code{Gender}).
#' @param visit A character string or factor specifying a single value in the
#'   \code{Visit} column to filter the data by. If \code{NULL} (default), no
#'   filtering is applied based on \code{Visit}.
#' @param y_variable A **bare column name** (non-standard evaluation) representing
#'   the numeric variable to be plotted on the y-axis (e.g., \code{Value}, \code{BMI}).
#' @param y_label A character string for the y-axis label.
#' @param plot_title A character string for the main title of the generated plot.
#' @param y_limits A numeric vector of length 2 (e.g., \code{c(min, max)}) to
#'   set the limits of the y-axis. Defaults to \code{NULL} (uses \code{ggplot2} defaults).
#' @param table A logical value. If \code{TRUE} (default), the function combines
#'   the box plot with a summary statistics table using \code{\link[ggpubr]{ggsummarystats}}.
#'   If \code{FALSE}, only the plot is returned.
#'
#' @return Returns a \code{ggplot} object or a list of plots/tables (a \code{ggpubr}
#'   object) if \code{table = TRUE}. The output displays the box plot, potentially
#'   combined with a summary table.
#'
#' @details
#' The function relies heavily on functions from the \strong{ggpubr} package,
#' specifically \code{\link[ggpubr]{ggboxplot}} and \code{\link[ggpubr]{ggsummarystats}}.
#'
#' \itemize{
#'   \item **Box Plot Aesthetics:** The plot uses the \code{"lancet"} color palette,
#'         includes individual data points using \code{add = "jitter"}, and overlays
#'         the **mean** as a large green point (\code{\#00AC10}). 
#'   \item **Non-Standard Evaluation (NSE):** The \code{group} and \code{y_variable}
#'         arguments must be supplied as **bare column names**.
#'   \item **Data Validation:** The function includes checks to ensure the specified
#'         \code{y_variable} column exists and is numeric, stopping execution if not.
#' }
#'
#' @importFrom dplyr filter %>%
#' @importFrom ggpubr ggboxplot stat_summary ggsummarystats
#' @importFrom ggplot2 labs theme_minimal theme element_text element_blank ylim
#'
#' @examples
#' \dontrun{
#' # Assuming 'baseline_data' has columns: 'RandArm', 'Age', and 'Visit'
#'
#' # 1. Plot Age (numeric) grouped by RandArm (categorical), showing plot + table
#' plot_and_table <- boxplot_with_table1(
#'   data = baseline_data,
#'   group = RandArm,
#'   y_variable = Age,
#'   y_label = "Age (Years)",
#'   plot_title = "Age Distribution at Baseline",
#'   visit = "Baseline"
#' )
#' print(plot_and_table)
#'
#' # 2. Plot only the box plot (no summary table)
#' only_plot <- boxplot_with_table1(
#'   data = baseline_data,
#'   group = RandArm,
#'   y_variable = Age,
#'   y_label = "Age (Years)",
#'   plot_title = "Age Distribution at Baseline",
#'   visit = "Baseline",
#'   table = FALSE
#' )
#' print(only_plot)
#' }
#'
boxplot_with_table1 <- function(data, group, visit = NULL, y_variable, y_label, plot_title, y_limits = NULL, table = TRUE) {
  
  # Convert variable names to strings
  group_name <- deparse(substitute(group))
  y_variable_name <- deparse(substitute(y_variable))
  
  # Filter dataset for a specific visit if needed
  if (!is.null(visit)) {
    data <- data %>% filter(Visit == visit)
  }
  
  # Check column existence
  if (!(y_variable_name %in% colnames(data))) {
    stop(paste("Error: Column", y_variable_name, "not found in dataset."))
  }
  
  # Ensure the y_variable column is numeric
  if (!is.numeric(data[[y_variable_name]])) {
    stop(paste("Error: Column", y_variable_name, "must be numeric."))
  }
  
  # Create the boxplot
  plot <- ggboxplot(
    data,
    x = group_name,
    y = y_variable_name,
    color = group_name,
    fill = group_name,
    palette = "lancet",
    alpha = 0.2,
    add = "jitter"
  ) +
    stat_summary(fun = mean, geom = "point", shape = 20, size = 3.5, color = "#00AC10") +
    labs(x = NULL, y = y_label, title = plot_title) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "none"
    )
  
  if (!is.null(y_limits)) {
    plot <- plot + ylim(y_limits)
  }
  
  # Generate summary table only if table = TRUE
  if (table) {
    result <- ggsummarystats(
      data = data,
      x = group_name,
      y = y_variable_name,
      summaries = c("n", "mean", "sd", "median", "iqr"), digits = 1,
      ggfunc = function(...) plot
    )
    return(result)  # Return the combined plot + table
  }
  
  return(plot)  # Return only the plot if table = FALSE
}