#' @title Plot ROC Curve, Optimal Threshold, and AUC CI
#'
#' @description
#' Calculates the Receiver Operating Characteristic (ROC) curve, determines the 
#' Area Under the Curve (AUC) with its 95\% Confidence Interval (CI) using 
#' DeLong's method, and finds the optimal classification threshold using Youden's 
#' J statistic. It then generates a highly customized plot visualizing these results.
#'
#' @details
#' This function requires the \code{pROC} and \code{tidyverse} packages (specifically \code{dplyr}). 
#' The X-axis of the plot represents the False Positive Rate (1 - Specificity), 
#' and the Y-axis represents the True Positive Rate (Sensitivity). 
#' The optimal threshold is marked and annotated with its corresponding \code{marker_at_predfrom} value.
#' The logic for determining the \code{marker_at_predfrom} value is to find the highest predicted risk 
#' score that is less than or equal to the optimal Youden's threshold.
#'
#' @param df A data frame containing at least two columns: the actual class labels 
#'   (must be named \code{status}) and the predicted risk scores. Must also contain 
#'   a column named \code{marker_at_predfrom} for annotation.
#' @param what The name of the column in \code{df} that contains the predicted risk 
#'   probabilities or scores. Must be a numeric column. Defaults to \code{"Risk_PC"}.
#'
#' @return The function is executed for its side effect (generating the plot). It 
#'   invisibly returns the ROC object.
#' 
#' @importFrom pROC roc ci coords plot
#' @importFrom dplyr %>% arrange filter slice_tail
#' @importFrom graphics points text abline
#' @examples
#' \dontrun{
#' # Assuming 'status' is the true class and 'Risk_PC' is the prediction
#' library(pROC)
#' library(tidyverse)
#' 
#' set.seed(42)
#' example_data <- data.frame(
#'   status = sample(c(0, 1), 100, replace = TRUE),
#'   Risk_PC = runif(100),
#'   marker_at_predfrom = rnorm(100, mean = 50, sd = 10) 
#' )
#' 
#' plot_auc(example_data)
#' 
#' # Example with a different column name
#' example_data$Prediction <- example_data$Risk_PC
#' plot_auc(example_data, what = "Prediction")
#' }
#'
#'

plot_auc <- function(df,what = "Risk_PC") {
  
  df <- as.data.frame(df)
  roc_obj <- roc(df$status,df[,what],algorithm=1)
  plot(roc_obj, 
       main = "ROC Curve",
       col = "#0072B2", # Custom color for the curve
       lwd = 2,
       auc.polygon = TRUE, 
       grid = c(0.1, 0.2), 
       max.auc.polygon = TRUE, 
       auc.polygon.col = "lightblue")
  
  # Calculate the 95% CI using DeLong's method
  ci_delong <- ci(roc_obj, method = "delong")
  
  optimal_coords <- coords(roc_obj, x = "best", best.method = "youden")
  # Add the optimal point to the plot
  
  df$dup <- df[,what]
  df1 <- df %>% arrange(dup) %>% filter(dup<=optimal_coords$threshold) %>% slice_tail(n=1)
  points(
    x = optimal_coords$specificity, # X-axis is 1 - Specificity (FPR)
    y = optimal_coords$sensitivity,     # Y-axis is Sensitivity (TPR)
    pch = 19,                            # Character for a star
    col = "red", 
    cex = 1.5                           # Size of the star
  )
  
  # Optional: Add text label next to the star
  text(
    x = optimal_coords$specificity, 
    y = optimal_coords$sensitivity, 
    labels = paste0("Optimal Threshold\n (Marker = ", round(df1$marker_at_predfrom, 1), ")"),
    pos = 3, # Position 4 is to the right
    col = "red",
    cex = 0.8
  )
  
  p1 <- paste0("AUC: ", round(ci_delong[2], 3), "\n")
  p2 <- paste0("95% CI: [", round(ci_delong[1], 3), ", ", round(ci_delong[3], 3), "]\n")
  
  text(
    x = 0.5, 
    y = 0.5, 
    labels = paste(p1,p2),
    pos = 1, # Position 4 is to the right
    col = "blue",
    cex = 1.1
  )
  
}