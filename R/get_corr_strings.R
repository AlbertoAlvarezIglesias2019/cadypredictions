#' Compute and Format Spearman Correlation Results for a Set of Biomarkers
#'
#' This function computes the Spearman's rank correlation coefficient (\eqn{\rho}),
#' 95% confidence interval (CI), and p-value between a single specified
#' \code{predictor} variable and a fixed set of five baseline biomarkers. The
#' results are then formatted into presentation-ready strings using the helper
#' function \code{format_corr()}.
#'
#' @param data A data frame containing the predictor variable and all biomarker
#'   columns specified in the function's internal list.
#' @param predictor A character string specifying the name of the column to be
#'   used as the predictor variable (e.g., "Age", "BMI", "Ejection_Fraction").
#'
#' @return A named list of five character strings. Each element in the list
#'   corresponds to one of the five biomarkers, and its value is the formatted
#'   correlation result string (e.g., "Spearman's Rho = ..., 95% CI (...), p-value ...").
#'
#' @details
#' \itemize{
#'   \item **Biomarker Set**: The function is hardcoded to use the following
#'     five baseline biomarkers: \code{"BNP_bl"}, \code{"NT_pro_BNP_bl"},
#'     \code{"CRP_bl"}, \code{"hsTnI_STAT_bl"}, and \code{"Galectin_3_bl"}.
#'   \item **Correlation Method**: It uses \code{\link[psych]{corr.test}} with
#'     \code{method = "spearman"} (Spearman's rank correlation). No p-value
#'     adjustment is applied (\code{adjust = "none"}).
#'   \item **Formatting**: It relies on the presence of the helper function
#'     \code{format_corr()} (not defined here, but assumed to be available) to
#'     process the output of \code{psych::corr.test} into the final strings.
#' }
#'
#' @seealso
#' The helper function \code{\link{format_corr}} for the detailed string formatting.
#'
#' @importFrom dplyr select all_of %>%
#' @importFrom psych corr.test
#'
#' @examples
#' \dontrun{
#' # Assuming 'clinical_data' is a data frame with columns like 'Age' and the biomarker set.
#'
#' # Compute correlations between 'Age' and the five biomarkers
#' age_corr_strings <- get_corr_strings(
#'   data = clinical_data,
#'   predictor = "Age"
#' )
#'
#' # The result is a list like:
#' # $BNP_bl
#' # [1] "Spearman's Rho = 0.350, 95% CI (0.200-0.500), p-value 0.005"
#' # $NT_pro_BNP_bl
#' # [1] "Spearman's Rho = 0.420, 95% CI (0.280-0.560), p-value < 0.001"
#' # ...
#' }
#'
get_corr_strings <- function(data, predictor) {
  
  biomks <- c("BNP_bl","NT_pro_BNP_bl","CRP_bl","hsTnI_STAT_bl","Galectin_3_bl")
  
  res <- psych::corr.test(
    data %>% select(all_of(predictor), all_of(biomks)),
    method = "spearman",
    adjust = "none"
  )
  
  # Format correlation strings for rows 1–5
  setNames(
    lapply(1:5, function(i) format_corr(res, i)),
    biomks
  )
}