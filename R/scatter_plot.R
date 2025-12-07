#' Create a Single Scatter Plot with Linear Regression Line
#'
#' Generates a simple scatter plot comparing two numeric variables (\code{xvar}
#' and \code{yvar}) and overlays a **linear regression line (\code{geom_smooth(method = "lm")})**
#' with its 95% confidence interval (CI).
#'
#' @param data A data frame containing the variables specified by \code{xvar} and \code{yvar}.
#' @param xvar A character string specifying the column name for the variable to be
#'   plotted on the x-axis.
#' @param yvar A character string specifying the column name for the variable to be
#'   plotted on the y-axis.
#' @param subtitle A character string for the subtitle of the plot (often used for
#'   the correlation statistics).
#' @param ylab A character string for the y-axis label. This is also used as the main title.
#' @param ylim A numeric vector of length 2 (e.g., \code{c(min, max)}) to set the
#'   limits of the y-axis. Defaults to \code{NULL} (uses \code{ggplot2} defaults).
#'
#' @return A \code{ggplot} object displaying the scatter plot and the linear model fit.
#'
#' @details
#' \itemize{
#'   \item **Regression Line**: The line of best fit is calculated using the
#'         \strong{linear model (\code{lm})} method, and the standard error ribbon
#'         is shown by default. 

#[Image of a scatter plot with a linear regression line and a shaded confidence interval]

#'   \item **Labeling**: The \code{xvar} string is used directly for the x-axis label.
#'         The \code{ylab} string is used for both the y-axis label and the main plot title.
#'   \item **Non-Standard Evaluation (NSE)**: Column names are passed as character
#'         strings and accessed using the \code{.data} pronoun.
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth labs ylim
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' # Assuming 'patient_data' has columns 'BMI' and 'Cholesterol'
#'
#' corr_text <- "Spearman's Rho = 0.55, p-value < 0.001"
#'
#' plot_cholesterol <- scatter_plot(
#'   data = patient_data,
#'   xvar = "BMI",
#'   yvar = "Cholesterol",
#'   subtitle = corr_text,
#'   ylab = "Total Cholesterol (mg/dL)",
#'   ylim = c(100, 300)
#' )
#' print(plot_cholesterol)
#' }
#'
scatter_plot <- function(data, xvar, yvar, subtitle, ylab, ylim = NULL) {
  
  p <- ggplot(data, aes(x = .data[[xvar]], y = .data[[yvar]])) +
    geom_point(alpha = 0.8) +
    geom_smooth(method = "lm") +
    labs(x = xvar, y = ylab, title = ylab, subtitle = subtitle)
  
  if (!is.null(ylim)) p <- p + ylim(ylim)
  p
}