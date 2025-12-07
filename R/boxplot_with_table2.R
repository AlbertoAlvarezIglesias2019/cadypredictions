#' Create Box Plot with Summary Table Aggregated by Participant ID
#'
#' Generates a box plot for a specified numeric variable (\code{y_variable}),
#' grouped by a categorical variable (\code{group}). This function is designed
#' for **repeated measures data**. While the box plot shows the distribution of
#' all raw data points, the accompanying summary statistics table is calculated
#' using the **mean value per participant** (\code{participant_id}).
#'
#' @param data A data frame containing the data. It must include columns for
#'   \code{group}, \code{participant_id}, and \code{y_variable}.
#' @param group A **bare column name** (non-standard evaluation) representing the
#'   categorical variable to be plotted on the x-axis and used for grouping
#'   (e.g., \code{RandArm}, \code{Treatment}).
#' @param participant_id A **bare column name** (non-standard evaluation) representing
#'   the unique identifier for each subject (e.g., \code{SubjectID}). This column
#'   is used to aggregate the data before calculating summary statistics.
#' @param visit A character string or factor specifying a single value in the
#'   \code{Visit} column to filter the data by. If \code{NULL} (default), no
#'   filtering is applied based on \code{Visit}.
#' @param y_variable A **bare column name** (non-standard evaluation) representing
#'   the numeric variable to be plotted on the y-axis (e.g., \code{Value}, \code{Score}).
#' @param y_label A character string for the y-axis label.
#' @param plot_title A character string for the main title of the generated plot.
#' @param y_limits A numeric vector of length 2 (e.g., \code{c(min, max)}) to
#'   set the limits of the y-axis. Defaults to \code{NULL}.
#' @param table A logical value. If \code{TRUE} (default), the function combines
#'   the box plot with a summary statistics table (n, mean, sd, median, IQR)
#'   derived from the aggregated participant means. If \code{FALSE}, only the plot is returned.
#'
#' @return Returns a \code{ggplot} object or a list of plots/tables (a \code{ggpubr}
#'   object) if \code{table = TRUE}.
#'
#' @details
#' This function is crucial for data where multiple measurements per subject exist
#' (e.g., multiple time points).
#'
#' \itemize{
#'   \item **Box Plot:** Displays the **raw distribution** of all measurements.
#'   \item **Summary Table:** Calculates a mean value for \code{y_variable} for
#'         each unique \code{participant_id} *before* calculating the overall group
#'         summary statistics (mean, median, SD, etc.). This prevents **pseudoreplication**
#'         in the summary table.
#'   \item The function uses \code{\link[ggpubr]{ggboxplot}} and
#'         \code{\link[ggpubr]{ggsummarystats}}.
#'   \item \code{group}, \code{participant_id}, and \code{y_variable} must be
#'         supplied as **bare column names**.
#' }
#'
#' @importFrom dplyr filter group_by summarise %>%
#' @importFrom ggpubr ggboxplot ggsummarystats
#' @importFrom ggplot2 labs theme_minimal theme element_text element_blank ylim
#' @importFrom rlang sym
#'
#' @examples
#' \dontrun{
#' # Assume 'longitudinal_data' has: 'RandArm', 'SubjectID', 'Value', 'Visit'
#'
#' # Plotting 'Value' grouped by 'RandArm' across all visits, with summaries
#' # calculated on the mean 'Value' per SubjectID.
#' plot_aggregated_summary <- boxplot_with_table2(
#'   data = longitudinal_data,
#'   group = RandArm,
#'   participant_id = SubjectID,
#'   y_variable = Value,
#'   y_label = "Observed Value (Units)",
#'   plot_title = "Value Distribution and Aggregated Summaries",
#'   table = TRUE
#' )
#' print(plot_aggregated_summary)
#' }
#'
boxplot_with_table2 <- function(data, group, participant_id, visit = NULL, y_variable, y_label, plot_title, y_limits = NULL, table = TRUE) {
  
  # Convert variable names to strings
  group_name <- deparse(substitute(group))
  y_variable_name <- deparse(substitute(y_variable))
  participant_id_name <- deparse(substitute(participant_id))
  
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
  
  # Aggregate by participant first
  data_summary <- data %>%
    group_by(!!sym(participant_id_name), !!sym(group_name)) %>%
    summarise(MeanValue = mean(!!sym(y_variable_name), na.rm = TRUE), .groups = "drop")
  
  # Boxplot using raw data
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
  
  # Return only the plot if table = FALSE
  if (!table) {
    return(plot)
  }
  
  # Summary table (based on participant-aggregated data)
  result <- ggsummarystats(
    data = data_summary,
    x = group_name,
    y = "MeanValue",
    summaries = c("n", "mean", "sd", "median", "iqr"), digits = 1,
    ggfunc = function(...) plot
  )
  
  return(result)  # Return the combined plot + table
}