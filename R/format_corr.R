#' Format Correlation Test Results into a String
#'
#' Takes the results from a correlation analysis (specifically, an object
#' structured like the output of \code{cor.test} or similar packages) and formats
#' the Spearman's rank correlation coefficient (\eqn{\rho}), its 95% confidence
#' interval (CI), and the p-value into a single, presentation-ready character string.
#'
#' @param res A list or object containing correlation results. It must contain a
#'   data frame named \code{ci} with columns for \code{r}, \code{lower}, \code{upper},
#'   and \code{p}.
#' @param row An integer specifying the row index within the \code{res$ci} data
#'   frame from which to extract the correlation statistics.
#' @param p_cutoff A numeric threshold. If the p-value is less than this cutoff,
#'   the p-value string will be formatted as \code{"< [p_cutoff]"}. Defaults to \code{0.001}.
#'
#' @return A character string containing the formatted correlation results,
#'   e.g., \code{"Spearman's Rho = 0.520, 95% CI (0.450-0.590), p-value < 0.001"}.
#'
#' @details
#' \itemize{
#'   \item **Rounding**: The correlation coefficient (\code{r}) and the CI limits
#'         are rounded to **three decimal places**. The p-value is also rounded
#'         to three decimal places unless it falls below the \code{p_cutoff}.
#'   \item **p-Value Handling**: If the calculated p-value is below the threshold
#'         (default is 0.001), it is displayed as \code{"< 0.001"} to prevent
#'         displaying excessive precision for very small values.
#' }
#'
#' @examples
#' \dontrun{
#' # Assuming 'my_corr_result' is an object containing correlation stats
#' # where my_corr_result$ci is a data frame
#'
#' # Example data structure for 'res'
#' res_example <- list(
#'   ci = data.frame(
#'     r = c(0.5197, 0.2301),
#'     lower = c(0.4503, 0.1509),
#'     upper = c(0.5891, 0.3102),
#'     p = c(0.000001, 0.015)
#'   )
#' )
#'
#' # Format the first row's results
#' formatted_string1 <- format_corr(res_example, row = 1)
#' # Result: "Spearman's Rho = 0.520, 95% CI (0.450-0.589), p-value < 0.001"
#'
#' # Format the second row's results
#' formatted_string2 <- format_corr(res_example, row = 2)
#' # Result: "Spearman's Rho = 0.230, 95% CI (0.151-0.310), p-value 0.015"
#' }
#'
format_corr <- function(res, row, p_cutoff = 0.001) {
  r  <- round(res$ci[row, "r"], 3)
  lo <- round(res$ci[row, "lower"], 3)
  up <- round(res$ci[row, "upper"], 3)
  p  <- res$ci[row, "p"]
  
  p_txt <- ifelse(p < p_cutoff, "< 0.001", round(p, 3))
  
  paste0("Spearman's Rho = ", r,
         ", 95% CI (", lo, "-", up, ")",
         ", p-value ", p_txt)
}