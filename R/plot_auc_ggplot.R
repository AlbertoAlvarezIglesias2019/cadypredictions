#' @title Plot Time-Dependent ROC Curve using ggplot2
#'
#' @description
#' Calculates the Receiver Operating Characteristic (ROC) curve and key performance 
#' metrics (AUC, CI, Optimal Threshold) and generates a publication-quality plot 
#' using the \code{ggplot2} framework. The plot visualizes the model's performance 
#' within a specified time window.
#'
#' @details
#' The function dynamically titles the plot based on the model type defined by the 
#' \code{what} argument (either "partly conditional" or "time dependent Cox PH"). 
#' The subtitle specifies the prediction time window using \code{predict_from} and 
#' \code{predict_to} columns.
#' 
#' The optimal threshold is determined by maximizing Youden's J statistic. The 
#' annotation for this point conditionally displays the corresponding 
#' \code{marker_at_predfrom} value if the marker and risk scores are monotonically 
#' related; otherwise, the raw probability threshold is used. The Area Under the 
#' Curve (AUC) polygon is correctly shaded, extending to the axes.
#'
#' @param df A data frame containing all necessary columns:
#'   \itemize{
#'     \item \code{status}: The binary actual outcome (0/1).
#'     \item \code{predict_from}, \code{predict_to}: Time points defining the prediction window (used in the subtitle).
#'     \item \code{marker_at_predfrom}: A numeric column representing the marker value at the time of prediction (used for threshold annotation).
#'     \item \code{marker_name}: The name of the biomarker being analyzed (used for plot annotation).
#'   }
#' @param what The name of the column in \code{df} that contains the predicted risk 
#'   probabilities or scores. Must be a numeric column. Defaults to \code{"Risk_PC"}.
#'
#' @return A \code{ggplot} object containing the ROC curve visualization.
#' 
#' @importFrom pROC roc ci coords
#' @importFrom dplyr %>% mutate arrange filter slice_tail bind_rows
#' @importFrom ggplot2 ggplot aes geom_polygon geom_line geom_abline geom_point 
#'             labs annotate coord_equal scale_x_continuous scale_y_continuous 
#'             theme_minimal theme element_text element_rect element_blank
#' @examples
#' \dontrun{
#' # Assuming necessary packages are loaded: library(pROC); library(tidyverse); library(ggplot2)
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
#' # Plotting for the "partly conditional" model
#' plot_auc_ggplot(example_data, what = "Risk_PC")
#' 
#' # Plotting for the "time dependent Cox PH" model
#' plot_auc_ggplot(example_data, what = "Risk_TD")
#' }
#'
#'
#'


plot_auc_ggplot <- function(df,what = "Risk_PC") {
  
  # Ensure input is a data frame
  df <- as.data.frame(df)
  
  # --- 1. Dynamic Title Creation ---
  if (what == "Risk_PC") {
    tt1 <- "ROC Curve using partly conditional models"
  } else {
    tt1 <- "ROC Curve using time dependent Cox PH models"
  }
  
  # Subtitle for the time window
  tt2 <- paste("CT prediction from", df$predict_from[1], "to", df$predict_to[1], "days")
  
  # --- 2. ROC Calculation (pROC) ---
  roc_obj <- roc(df$status, df[, what], algorithm = 1, quiet = TRUE)
  
  # Get coordinates for plotting
  roc_data <- coords(roc_obj, "all", ret = c("specificity", "sensitivity"))
  # Convert specificity to False Positive Rate (FPR) for ggplot X-axis
  roc_data <- roc_data %>% 
    mutate(FPR = 1 - specificity)
  
  # Calculate AUC and CI
  ci_delong <- ci(roc_obj, method = "delong")
  auc_text <- paste0("AUC = ", round(ci_delong[2], 3), "; 95% CI [", 
                     round(ci_delong[1], 3), ", ", round(ci_delong[3], 3), "]")
  
  # --- 3. Optimal Threshold Calculation (Youden's J) ---
  optimal_coords <- coords(roc_obj, x = "best", best.method = "youden")
  
  # Conditional Optimal Threshold Label Logic
  df$dup <- df[, what]
  temp1 <- order(df$marker_at_predfrom)
  temp2 <- order(df$dup)
  
  if (all(temp1 == temp2)) {
    # If ordering is the same, use the marker value
    df1 <- df %>% 
      arrange(dup) %>% 
      filter(dup <= optimal_coords$threshold) %>% 
      slice_tail(n = 1)
    threshold_label <- paste0("Optimal Threshold\n(Marker = ", round(df1$marker_at_predfrom, 1), ")")
  } else {
    # If ordering is different, use the raw probability threshold
    threshold_label <- paste0("Optimal Threshold\n(Prob = ", round(optimal_coords$threshold, 3), ")")
  }
  
  # Prepare optimal point data for plotting
  optimal_point_data <- data.frame(
    FPR = 1 - optimal_coords$specificity[1],
    Sensitivity = optimal_coords$sensitivity[1]
  )
  
  # Add the final closing point (1, 0) to the data.
  # This point is essential to close the polygon to the x-axis.
  closing_point <- data.frame(specificity = 0, sensitivity = 0, FPR = 1) 
  
  roc_data <- roc_data %>% arrange(FPR,sensitivity)
  
  polygon_data <- bind_rows(roc_data, closing_point) 
  
  # --- 4. ggplot Construction ---
  p <- ggplot(roc_data, aes(x = FPR, y = sensitivity)) +
    
    # 4a. Add the AUC Polygon (shaded area)
    geom_polygon(data = polygon_data, 
                 aes(x = FPR, y = sensitivity), 
                 fill = "lightblue", 
                 alpha = 0.5)  +
    
    # 4b. Add the ROC Curve Line
    geom_line(color = "#0072B2", linewidth = 0.7) +
    
    # 4c. Add the Diagonal Reference Line (Random Classifier)
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
    
    # 4d. Add Optimal Point Marker
    geom_point(data = optimal_point_data, aes(x = FPR, y = Sensitivity), 
               color = "red", size = 3, shape = 19) +
    
    # 4e. Add Labels and Annotations
    labs(
      title = tt1,
      subtitle = tt2,
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    
    # Annotate AUC and CI (Positioned near x=0.5, y=0.1)
    annotate("text", x = 0.5, y = 0.1, label = auc_text, color = "blue", size = 4.5) +
    
    # Annotate Optimal Threshold (Positioned above the point marker)
    annotate("text", x = optimal_point_data$FPR, y = optimal_point_data$Sensitivity, 
             label = threshold_label, color = "red", size = 3.5, vjust = -0.3) +
    
    # Annotate Marker Name (Positioned near top right corner)
    annotate("text", x = 0, y = 0.95, label = df$marker_name[1], color = "darkgreen", 
             size = 4.5,hjust=0) +
    
    # 4f. Set X and Y Axis Limits (0 to 1) and Aspect Ratio
    coord_equal() + # Ensures X and Y axes have the same scaling
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    
    # 4g. Apply a Clean Theme
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
    )
  
  return(p)
}