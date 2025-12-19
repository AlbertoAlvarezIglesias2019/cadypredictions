#' Log2 Transformation with Special Handling for Zeros
#'
#' This function performs a base-2 logarithmic transformation on a numeric vector. 
#' It handles zero values by substituting them with half of the minimum non-zero 
#' value found in the vector before applying the log transformation.
#'
#' @param x A numeric vector to be transformed.
#' @param na.rm Logical. Should missing values (including NaN) be removed when 
#' calculating the minimum non-zero value? Defaults to \code{TRUE} within the \code{min} call.
#'
#' @details 
#' In biological assays (like BNP or Troponin), a value of zero often indicates 
#' a result below the Limit of Detection (LOD). Since $log_2(0)$ is undefined ($-\infty$), 
#' this function imputes zero values using the "half-minimum" rule:
#' \enumerate{
#'   \item Identify the smallest value in the vector that is strictly greater than zero.
#'   \item Divide this value by 2 (the "replacement value").
#'   \item Replace all zeros in the original vector with this replacement value.
#'   \item Apply \code{log2()} to the adjusted vector.
#' }
#'
#' @return A numeric vector of the same length as \code{x}, containing the $log_2$ 
#' transformed values.
#'
#' @examples
#' marker_values <- c(0, 10, 20, 100, 0)
#' calc_log2_special(marker_values)
#' # Minimum non-zero is 10; zeros are replaced by 5 before log2 is applied.
#'
#' @export
#' 

calc_log2_special <- function(x) {
  # Find the minimum of the non-zero values
  # (Assuming you want the min of the whole column to define the 'floor')
  min_val <- min(x[x > 0], na.rm = TRUE)
  replacement_val <- min_val / 2
  
  # Apply the logic: if 0, use replacement_val; else use x
  adjusted_x <- ifelse(x == 0, replacement_val, x)
  
  return(log2(adjusted_x))
}