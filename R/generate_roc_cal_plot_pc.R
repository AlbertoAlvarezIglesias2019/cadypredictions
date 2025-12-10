#' @title Generate ROC Curve, Calibration Plot, and Metrics for a Binary Classifier (PC/TD Model)
#'
#' @description This function evaluates a binary classification model's discrimination
#' (using ROC curve and AUC) and calibration (using a calibration plot) on a
#' given dataset, typically derived from a partially conditional or time-dependent (TD)
#' model. It also calculates an optimal cutpoint, the resulting confusion matrix,
#' and a Wilson confidence interval for the proportion of correctly diagnosed positive cases.
#'
#' @details The function performs the following steps:
#' 1. **Optimal Cutpoint:** Calculates the optimal cutpoint using the \code{cutpointr} package based on the specified method and metric (default: maximize the sum of sensitivity and specificity).
#' 2. **ROC Curve:** Generates the ROC curve using \code{plotROC} and annotates it with the AUC value and the point corresponding to the optimal cutpoint. 
#' 3. **Calibration Plot:** Creates a calibration plot by grouping predicted risks into deciles and plotting the mean predicted risk vs. the mean observed event rate in each group. The plot includes a line of perfect calibration ($y=x$). 
#' 4. **Plot Merging:** Combines the ROC and Calibration plots side-by-side using \code{patchwork}.
#' 5. **Metric Calculation:** Calculates the confusion matrix using \code{caret} based on the optimal cutpoint, and calculates a one-sided 95% Wilson confidence interval for the true positive rate (sensitivity).
#'
#' @param data A data frame containing the risk predictions and the observed event status.
#' @param risk_col A string specifying the name of the column in \code{data} that contains the numeric probability predictions (risk scores). Default is \code{"Risk_PC"}.
#' @param event_col A string specifying the name of the column in \code{data} that contains the observed binary event status (0 or 1). Default is \code{"status"}.
#' @param cutpointr_method A function specifying the method to find the optimal cutpoint, as used by \code{cutpointr::cutpointr}. Default is \code{cutpointr::maximize_metric}.
#' @param cutpointr_metric A function specifying the metric to maximize when finding the optimal cutpoint, as used by \code{cutpointr::cutpointr}. Default is \code{cutpointr::sum_sens_spec}.
#' @param cutpointr_boots An integer specifying the number of bootstrap runs for \code{cutpointr::cutpointr} to estimate variability. Default is \code{1000}.
#' @param roc_title A string for the title of the ROC plot. Default is \code{"ROC Curve"}.
#' @param cal_title A string for the title of the calibration plot. Default is \code{"Calibration Plot"}.
#' @param merged_title A string for the title of the merged plot. Default is \code{"Partly conditional model discrimination and calibration (on test set)"}.
#'
#' @return A list with the following components:
#' \item{AUC_Value}{The Area Under the ROC Curve (AUC).}
#' \item{Optimal_Cutpoint}{The calculated optimal cutpoint.}
#' \item{ROC_Curve}{A \code{ggplot} object of the ROC curve.}
#' \item{Calibration_Plot}{A \code{ggplot} object of the calibration plot.}
#' \item{Merged_Plot}{A \code{patchwork} object combining the ROC and Calibration plots.}
#' \item{Confusion_Matrix}{The \code{confusionMatrix} object from the \code{caret} package using the optimal cutpoint.}
#' \item{Wilson_CI_Output}{A formatted string containing the mean proportion correctly diagnosed (sensitivity/True Positive Rate) and its one-sided 95% Wilson confidence interval.}
#'
#' @importFrom rlang sym
#' @importFrom dplyr mutate group_by summarise n ungroup
#' @importFrom stats quantile binom.confint
#' @import cutpointr
#' @import plotROC
#' @import ggplot2
#' @import patchwork
#' @import caret
#'
#' @examples
#' \dontrun{
#' # Assuming 'my_data' has 'Risk_PC' and 'status' columns:
#' # my_results <- generate_roc_cal_plot_pc_td(data = my_data)
#' # print(my_results$Merged_Plot)
#' # print(my_results$Confusion_Matrix)
#' }
#'
generate_roc_cal_plot_pc_td <- function(
    data,
    risk_col = "Risk_PC",
    event_col = "status",
    cutpointr_method = maximize_metric,
    cutpointr_metric = sum_sens_spec,
    cutpointr_boots = 1000,
    roc_title = "ROC Curve",
    cal_title = "Calibration Plot",
    merged_title = "Partly conditional model discrimination and calibration (on test set)"
) {
  
  # --------------------------------------------------------------------
  # 1. INPUT DATA PREPARATION
  # --------------------------------------------------------------------
  out <- data
  
  # --------------------------------------------------------------------
  # 2. OPTIMAL CUTPOINT (cutpointr)
  # --------------------------------------------------------------------
  cp <- cutpointr::cutpointr(
    data = out,
    x = !!rlang::sym(risk_col),      # numeric probability predictions
    class = !!rlang::sym(event_col), # FIXED tidy-eval string column name
    method = cutpointr_method,
    metric = cutpointr_metric,
    pos_class = 1,
    neg_class = 0,
    direction = ">=",
    boot_runs = cutpointr_boots
  )
  
  optimal_cutpoint <- cp$optimal_cutpoint[1]
  
  # --------------------------------------------------------------------
  # 3. ROC CURVE (using plotROC)
  # --------------------------------------------------------------------
  basicplot <- ggplot2::ggplot(
    out,
    aes(
      d = !!rlang::sym(event_col),
      m = !!rlang::sym(risk_col)
    )
  ) +
    plotROC::geom_roc()
  
  auc_val <- plotROC::calc_auc(basicplot)$AUC
  
  roc_plot <- basicplot +
    plotROC::style_roc(theme = ggplot2::theme_grey) +
    ggplot2::theme(axis.text = ggplot2::element_text(colour = "blue")) +
    ggplot2::ggtitle(roc_title) +
    ggplot2::annotate("text",
                      x = .75, y = .25,
                      label = paste("AUC =", round(auc_val, 3))) +
    ggplot2::scale_x_continuous("1 - Specificity",
                                breaks = seq(0, 1, 0.1)) +
    # Add the optimal cutpoint point in ROC space
    ggplot2::geom_point(aes(x = 1 - cp$specificity,
                            y = cp$sensitivity),
                        color = "red",
                        size = 3) +
    ggplot2::annotate("text",
                      x = .75, y = .20,
                      label = paste("Optimal cutpoint =", round(optimal_cutpoint, 4)))
  
  # --------------------------------------------------------------------
  # 4. CALIBRATION PLOT
  # --------------------------------------------------------------------
  out_cal <- out %>%
    dplyr::mutate(
      Groups = cut(
        !!rlang::sym(risk_col),
        stats::quantile(out[[risk_col]], prob = seq(0, 1, 0.1)),
        include.lowest = TRUE
      )
    )
  
  tmpDat <- out_cal %>%
    dplyr::group_by(Groups) %>%
    dplyr::summarise(
      X = mean(.data[[risk_col]]),
      Y = mean(.data[[event_col]]),
      N = dplyr::n()
    ) %>%
    dplyr::ungroup()
  
  cal_plot <- ggplot2::ggplot(tmpDat, aes(x = X, y = Y)) +
    ggplot2::geom_abline(intercept = 0, slope = 1,
                         color = "red", size = 2) +
    ggplot2::geom_point(shape = 1, stroke = 3) +
    ggplot2::scale_size_continuous(range = c(2, 20)) +
    ggplot2::ggtitle(cal_title) +
    ggplot2::xlab("Predicted event rate") +
    ggplot2::ylab("Observed event rate") +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 14),
      axis.title = ggplot2::element_text(size = 14, face = "bold")
    )
  
  # --------------------------------------------------------------------
  # 5. MERGE PLOTS (patchwork)
  # --------------------------------------------------------------------
  merged_plot <- roc_plot + cal_plot +
    patchwork::plot_annotation(title = merged_title)
  
  # --------------------------------------------------------------------
  # 6. CONFUSION MATRIX + WILSON CI
  # --------------------------------------------------------------------
  cm_data <- out %>%
    dplyr::mutate(
      pred_class_prob = if_else(.data[[risk_col]] >= optimal_cutpoint, 1, 0),
      pred_class_prob = factor(pred_class_prob, levels = c(0, 1)),
      observed = factor(.data[[event_col]], levels = c(0, 1))
    )
  
  cm_prob <- caret::confusionMatrix(
    data = cm_data$pred_class_prob,
    reference = cm_data$observed,
    positive = "1"
  )
  
  cm_tab <- cm_prob$table
  TP <- cm_tab["1", "1"]
  FN <- cm_tab["0", "1"]
  n_CT <- TP + FN
  
  wilson_res <- binom.confint(
    x = TP, n = n_CT,
    methods = "wilson",
    conf.level = 0.90     # 90% → equivalent to 95% one-sided lower
  )
  
  wilson_text <- paste0(
    "Mean proportion correctly diagnosed = ", round(wilson_res$mean, 2), ", ",
    "one-sided 95% CI (", round(wilson_res$lower, 2), " - ",
    round(wilson_res$upper, 2), ")"
  )
  
  # --------------------------------------------------------------------
  # 7. RETURN RESULTS
  # --------------------------------------------------------------------
  return(list(
    AUC_Value = auc_val,
    Optimal_Cutpoint = optimal_cutpoint,
    ROC_Curve = roc_plot,
    Calibration_Plot = cal_plot,
    Merged_Plot = merged_plot,
    Confusion_Matrix = cm_prob,
    Wilson_CI_Output = wilson_text
  ))
}