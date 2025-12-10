This function creates a classic calibration plot (also known as a calibration curve or reliability plot) using ggplot2. Here is the Roxygen documentation for your plot_calibration_ggplot function.
📝 Roxygen Documentation for plot_calibration_ggplot
R

#' @title Plot Model Calibration Curve using ggplot2
#'
#' @description
#' Calculates and plots a calibration curve for predicted event rates against observed 
#' event rates within a specified time window, using the \code{ggplot2} framework. 
#' This plot is crucial for assessing the agreement between predicted probabilities 
#' and actual outcomes (calibration).
#'
#' @details
#' The calibration process involves the following steps:
#' \enumerate{
#'   \item The predicted risk scores (\code{what} column) are grouped into 10 bins 
#'         using quantiles (\code{seq(0, 1, 0.1)}).
#'   \item The observed outcome (\code{observed}) is determined by checking if the 
#'         event occurred (\code{status}) within the prediction window 
#'         (\code{predict_from} to \code{predict_to}).
#'   \item For each group, the mean predicted risk (\code{X}) and the mean observed 
#'         rate (\code{Y}) are calculated.
#' }
#' The plot displays these binned points against the ideal identity line (perfect 
#' calibration, slope=1). The closer the points are to the line, the better the 
#' model's calibration. The plot titles dynamically reflect the model type and the 
#' prediction time window. 
#'
#' @param df A data frame containing all necessary columns:
#'   \itemize{
#'     \item \code{status}: The binary actual outcome (0/1).
#'     \item \code{time_to_event}: The time until the event occurred (or censoring time).
#'     \item \code{predict_from}, \code{predict_to}: Time points defining the prediction window (used for defining the 'observed' outcome and in the subtitle).
#'   }
#' @param what The name of the column in \code{df} that contains the predicted risk 
#'   probabilities or scores. Must be a numeric column. Defaults to \code{"Risk_PC"}.
#'
#' @return A \code{ggplot} object containing the calibration curve visualization.
#' 
#' @importFrom dplyr %>% as_data_frame mutate group_by summarise ungroup if_else
#' @importFrom ggplot2 ggplot aes geom_abline geom_point scale_size_continuous 
#'             labs theme element_text xlim ylim
#' @examples
#' \dontrun{
#' # Assuming necessary packages are loaded: library(tidyverse); library(ggplot2)
#' set.seed(42)
#' example_data <- data.frame(
#'   status = sample(c(0, 1), 100, replace = TRUE, prob=c(0.7, 0.3)),
#'   time_to_event = runif(100, 100, 2000),
#'   Risk_PC = runif(100, 0, 0.5), # Predictions (must be somewhat close to observed rates)
#'   predict_from = rep(365, 100),
#'   predict_to = rep(1095, 100)
#' )
#' 
#' # Plotting for the "partly conditional" model
#' plot_calibration_ggplot(example_data, what = "Risk_PC")
#' 
#' # Example with a different column name
#' example_data$Risk_TD <- example_data$Risk_PC * 1.1
#' plot_calibration_ggplot(example_data, what = "Risk_TD")
#' }
#'
#'

plot_calibration_ggplot <- function(df,what = "Risk_PC") {
  
  # Ensure input is a data frame
  df <- as.data.frame(df)
  
  # --- 1. Dynamic Title Creation ---
  if (what == "Risk_PC") {
    tt1 <- "Calibration using partly conditional models (test set)"
  } else {
    tt1 <- "Calibration using time dependent Cox PH models (test set)"
  }
  
  # Subtitle for the time window
  tt2 <- paste("CT prediction from", df$predict_from[1], "to", df$predict_to[1], "days")
  
  df$dup <- df[, what]
  
  
  ################
  ### Calibration
  ################
  df <- df %>%
    mutate(Groups = cut(dup,quantile(df$dup,prob = seq(0,1,0.1)),include.lowest=TRUE)) %>% 
    mutate(observed = if_else(time_to_event>predict_to | time_to_event<predict_from,0,status))
  
  tmpDat <- df %>%
    group_by(Groups) %>%
    summarise(X = mean(dup),
              Y = mean(observed),
              N = n()) %>%
    ungroup()
  
  
  ####################
  ### Create the plot
  ####################
  xm <- max(tmpDat$X)
  ym <- max(tmpDat$Y)
  ggplot(as.data.frame(tmpDat), aes(x=X, y=Y)) +
    geom_abline(intercept = 0, slope = 1,color="red",size=2)+
    geom_point(shape=1,stroke=3) +     # Use hollow circles
    scale_size_continuous(range = c(2, 20)) +
    xlab("Predicted event rate") +
    ylab("Observed event rate") + 
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=14,face="bold"),
          legend.title = element_text(size=20),
          legend.text=element_text(size=20),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5)) +
    xlim(0,xm) + ylim(0,ym) + 
    labs(
      title = tt1,
      subtitle = tt2)
  
  
}