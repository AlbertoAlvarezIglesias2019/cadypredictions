#' @title Generate ROC, Calibration Plots, and Diagnostic Metrics for Survival Model
#'
#' @description
#' Calculates risk scores from a pre-fitted Cox model, performs cutpoint optimization,
#' generates an ROC curve (using plotROC), a calibration plot, and computes a
#' Confusion Matrix and diagnostic accuracy metrics.
#'
#' @param data A data frame (e.g., the test set) for which the risk scores should be calculated.
#' @param coxph_model An object of class \code{\link[survival]{coxph}} that has already been fitted.
#' @param time_col A character string representing the name of the time-to-event column in \code{data} (e.g., "time_to_event").
#' @param event_col A character string representing the name of the event/status column in \code{data} (e.g., "status").
#' @param survC_type A character string for the 'type' argument in \code{survC::calc_risk_score()}. Default is "risk".
#' @param tdroc_t The required time point 't' (numeric) for AUC calculation using \code{survC::tdroc_calc()}.
#' @param cutpointr_method The method used to select the optimal cutpoint in \code{\link[cutpointr]{cutpointr}}.
#' @param cutpointr_metric The metric to be optimized in \code{\link[cutpointr]{cutpointr}}.
#' @param cutpointr_boots The number of bootstrap runs for \code{\link[cutpointr]{cutpointr}}.
#' @param roc_title A character string for the title of the ROC curve plot.
#' @param cal_title A character string for the title of the calibration plot.
#' @param merged_title A character string for the overall title of the merged plot.
#'
#' @return A list containing the requested outputs:
#' \item{AUC_Value}{The AUC value obtained from \code{survC::tdroc_calc()} at time \code{tdroc_t}.}
#' \item{Optimal_Cutpoint}{The optimal risk score cutpoint determined by \code{cutpointr}.}
#' \item{ROC_Curve}{The ggplot object for the ROC curve with cutpoint highlighted.}
#' \item{Calibration_Plot}{The ggplot object for the calibration plot.}
#' \item{Merged_Plot}{The combined ROC and Calibration plot (using \code{patchwork/cowplot}).}
#' \item{Confusion_Matrix}{The output of \code{caret::confusionMatrix}.}
#' \item{Wilson_CI_Output}{A formatted character string with the Wilson CI for diagnostic accuracy.}
#'
generate_roc_cal_plot <- function(
    data,
    coxph_model, 
    time_col,    
    event_col,   
    survC_type = "risk", 
    tdroc_t,             
    cutpointr_method = maximize_metric, 
    cutpointr_metric = sum_sens_spec,   # Changed default to match raw code
    cutpointr_boots = 1000,             
    roc_title = "ROC Curve",            # Changed default to match raw code
    cal_title = "Calibration plot",     # Changed default to match raw code
    merged_title = "Cox PH model discrimination and calibration (on test set)"
) {
  
  # Make a mutable copy of the data
  data_out <- data
  
  # --- 1. Predictions ---
  # Calculate risk score (numeric vector) and add to data
  data_out$risk_coxph <- survC::calc_risk_score(
    model = coxph_model,
    data = data_out,
    type = survC_type
  )
  
  # --- 2. Data Filtering and AUC Calculation (tdroc_calc) ---
  # Filter data for non-missing values
  out <- data_out %>%
    filter(!is.na(.data[[time_col]]),
           !is.na(.data[[event_col]]),
           !is.na(risk_coxph))
  
  # AUC Calculation
  auc_240_res <- survC::tdroc_calc(
    time   = out[[time_col]],
    status = out[[event_col]],
    marker = out$risk_coxph,
    t      = tdroc_t
  )
  # Extract the AUC value (assumes tdroc_calc returns a list/vector where AUC is at index 2)
  auc_value <- auc_240_res$AUC[2] 
  
  # --- 3. Optimal Cutpoint (cutpointr) ---
  cp <- cutpointr::cutpointr(
    data = out,
    x = risk_coxph,
    class = .data[[event_col]], # Use event_col to reference the status column
    method = cutpointr_method,
    metric = cutpointr_metric,
    pos_class = 1,
    neg_class = 0,
    direction = ">=",
    boot_runs = cutpointr_boots
  )
  
  optimal_cutpoint <- cp$optimal_cutpoint[1]
  
  # --- 4. ROC Curve Plot (using plotROC) ---
  basicplot <- ggplot2::ggplot(out , aes(d = .data[[event_col]], m = risk_coxph)) + 
    plotROC::geom_roc()
  
  # NOTE: The raw code uses patchwork/cowplot style for merging, assuming availability.
  roc_plot <- basicplot + 
    plotROC::style_roc(theme = ggplot2::theme_grey) +
    ggplot2::theme(axis.text = ggplot2::element_text(colour = "blue")) +
    ggplot2::ggtitle(roc_title) + 
    ggplot2::annotate("text", 
                      x = .75, 
                      y = .25, 
                      label = paste("AUC =", round(auc_value, 3))) +
    ggplot2::scale_x_continuous("1 - Specificity", breaks = seq(0, 1, by = .1)) + 
    ggplot2::geom_point(aes(x = 1 - cp$specificity,
                            y = cp$sensitivity),
                        color = "red",
                        size = 3) +
    ggplot2::annotate("text",
                      x = .75,
                      y = .20,
                      label = paste("Optimal cutpoint =", round(optimal_cutpoint, 4)))
  
  # --- 5. Calibration Plot ---
  # Create the groups (must use the filtered 'out' dataset)
  out <- out %>%
    dplyr::mutate(
      Groups = cut(risk_coxph, 
                   stats::quantile(out$risk_coxph, prob = seq(0, 1, 0.1)), 
                   include.lowest = TRUE))
  
  tmpDat <- out %>%
    dplyr::group_by(Groups) %>%
    dplyr::summarise(X = mean(risk_coxph),
                     Y = mean(.data[[event_col]]), # Use event_col for status mean
                     N = dplyr::n()) %>%
    dplyr::ungroup()
  
  # Create the plot
  cal_plot <- ggplot2::ggplot(as.data.frame(tmpDat), aes(x = X, y = Y)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "red", size = 2) +
    ggplot2::geom_point(shape = 1, stroke = 3) + 
    ggplot2::scale_size_continuous(range = c(2, 20)) +
    ggplot2::ggtitle(cal_title) + 
    ggplot2::xlab("Predicted event rate") +
    ggplot2::ylab("Observed event rate") + 
    ggplot2::theme(axis.text = ggplot2::element_text(size = 14),
                   axis.title = ggplot2::element_text(size = 14, face = "bold"),
                   legend.title = ggplot2::element_text(size = 20),
                   legend.text = ggplot2::element_text(size = 20)) 
  
  # --- 6. Merging ROC curve and calibration plots ---
  # Uses the '+' operator which typically implies the 'patchwork' package, 
  # but cowplot::plot_grid is safer if patchwork is not available.
  # Assuming 'patchwork' for the '+' operator as per the raw code's final line.
  merged_plot <- roc_plot + cal_plot + 
    patchwork::plot_annotation(title = merged_title)
  
  # --- 7. Confusion matrix and Wilson CI ---
  opt_cp <- optimal_cutpoint # Already extracted
  
  cm_data <- out %>% 
    dplyr::mutate(
      pred_class_prob = if_else(risk_coxph >= opt_cp, 1, 0),
      pred_class_prob = factor(pred_class_prob, levels = c(0, 1)),
      observed = factor(.data[[event_col]], levels = c(0, 1))
    )
  
  conf_matrix <- caret::confusionMatrix(
    data = cm_data$pred_class_prob,
    reference = cm_data$observed,
    positive = "1"
  )
  
  # Percentage of primary biomarker diagnoses that match the definition of cardiotoxicity
  cm_tab <- conf_matrix$table
  
  TP <- cm_tab["1","1"]
  FN <- cm_tab["0","1"]
  
  n_CT <- TP + FN
  
  # Wilson's confidence interval
  wilson_res <- binom.confint(
    x = TP,
    n = n_CT,
    methods = "wilson",
    conf.level = 0.90
  )
  
  wilson_output_string <- paste0(
    "Mean proportion correctly diagnosed = ", round(wilson_res$mean, 2), ", ",
    "one-sided 95% CI (", round(wilson_res$lower, 2), " - ",  
    round(wilson_res$upper, 2), ")"
  )
  
  # --- 8. Prepare Output ---
  return(list(
    AUC_Value = auc_value,
    Optimal_Cutpoint = optimal_cutpoint,
    ROC_Curve = roc_plot,
    Calibration_Plot = cal_plot,
    Merged_Plot = merged_plot,
    Confusion_Matrix = conf_matrix,
    Wilson_CI_Output = wilson_output_string
  ))
}