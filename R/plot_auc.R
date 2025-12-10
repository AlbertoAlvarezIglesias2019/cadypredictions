#' @title Plot Time-Dependent ROC Curve and Optimal Threshold
#'
#' @description
#' Calculates and plots the Receiver Operating Characteristic (ROC) curve for 
#' a binary classification model (often from survival analysis), determines the 
#' Area Under the Curve (AUC) with its 95\% Confidence Interval (CI) using 
#' DeLong's method, and finds the optimal classification threshold using Youden's 
#' J statistic.
#'
#' @details
#' The main title of the plot is dynamically generated based on the \code{what} 
#' argument, distinguishing between models using "partly conditional" logic 
#' (\code{"Risk_PC"}) and "time-dependent Cox PH" logic (any other value). 
#' The title also includes the prediction time window defined by the \code{predict_from} 
#' and \code{predict_to} columns.
#' 
#' The optimal threshold is marked on the plot. The label for this threshold 
#' conditionally displays either the corresponding \code{marker_at_predfrom} value 
#' or the raw probability threshold. This condition checks if the ordering of the 
#' risk scores (\code{what} column) is the same as the ordering of the 
#' \code{marker_at_predfrom} values. If the orderings match, the marker value is used 
#' for annotation; otherwise, the raw probability is used.
#'
#' @param df A data frame containing:
#'   \itemize{
#'     \item \code{status}: The binary actual outcome (0/1).
#'     \item \code{predict_from}, \code{predict_to}: Time points defining the prediction window (used in the subtitle).
#'     \item \code{marker_at_predfrom}: A numeric column representing the marker value at the time of prediction (used for annotation).
#'     \item \code{marker_name}: The name of the biomarker being analyzed (used for annotation).
#'   }
#' @param what The name of the column in \code{df} that contains the predicted risk 
#'   probabilities or scores. Must be a numeric column. Defaults to \code{"Risk_PC"}.
#'
#' @return The function is executed for its side effect (generating the plot). It 
#'   invisibly returns the ROC object.
#' 
#' @importFrom pROC roc ci coords plot
#' @importFrom dplyr %>% arrange filter slice_tail
#' @importFrom graphics points text
#' @examples
#' \dontrun{
#' # Assuming necessary packages are loaded: library(pROC); library(tidyverse)
#' set.seed(42)
#' example_data <- data.frame(
#'   status = sample(c(0, 1), 100, replace = TRUE),
#'   Risk_PC = runif(100),
#'   Risk_TD = runif(100),
#'   predict_from = rep(365, 100),
#'   predict_to = rep(1095, 100),
#'   marker_at_predfrom = rnorm(100, mean = 50, sd = 10),
#'   marker_name = rep("C-Reactive Protein", 100)
#' )
#' 
#' # 1. Plotting for the "partly conditional" model
#' plot_auc(example_data, what = "Risk_PC")
#' 
#' # 2. Plotting for the "time dependent Cox PH" model
#' plot_auc(example_data, what = "Risk_TD")
#' }
#'
#'
#'

plot_auc <- function(df,what = "Risk_PC") {
  
  if (what == "Risk_PC") {
    tt1 <- paste("ROC Curve using partly contidional models")
  } else {
    tt1 <- paste("ROC Curve using time dependent Cox PH models")
  }
  
  tt2 <- paste("(CT prediction from ",df$predict_from[1]," to ",df$predict_to[1]," days)",sep="")
  
  tt <- paste(tt1,"\n",tt2,sep="")
  
  df <- as.data.frame(df)
  roc_obj <- roc(df$status,df[,what],algorithm=1)
  plot(roc_obj,
       main = tt,
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
  
  temp1 <- order(df$marker_at_predfrom)
  temp2 <- order(df$dup)
  
  
  if (all(temp1==temp2)) {
    df1 <- df %>% arrange(dup) %>% filter(dup<=optimal_coords$threshold) %>% slice_tail(n=1)
    tex <- paste0("Optimal Threshold\n(Marker = ", round(df1$marker_at_predfrom, 1), ")")
  } else {
    tex <- paste0("Optimal Threshold\n(Prob = ", round(optimal_coords$threshold, 3), ")")
  }
  
  
  
  
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
    labels = tex,
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
  
  text(
    x = 1, 
    y = 0.95, 
    labels = df$marker_name[1],
    pos = 4, # Position 4 is to the right
    col = "darkgreen",
    cex = 1.1
  )
  
}
