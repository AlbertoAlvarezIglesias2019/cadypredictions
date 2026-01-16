#' Calculate AUC and Classification Metrics with Wilson CIs
#'
#' @description
#' This function systematically evaluates the predictive performance of biomarker 
#' combinations across multiple clinical endpoints. It fits six survival models per 
#' combination and calculates the Area Under the Curve (AUC) and diagnostic 
#' metrics (Sensitivity, Specificity, PPV, NPV) at the optimal Youden's J threshold. 
#' Significantly, it incorporates Wilson Score intervals for sensitivity and 
#' DeLong confidence intervals for AUC.
#'
#' @param dir_in Character. Path to input CSVs (endpoints and \code{baseline_data.csv}).
#' @param dir_out Character. Path for saving \code{.RData} and \code{.csv} results.
#' @param dat_nam Character vector. Base names of endpoint data files.
#' @param mar_nam Character vector. Pool of biomarkers to generate combinations from.
#' @param predictores Character vector. Covariates for Adjusted models.
#' @param saveopt Character (optional). Prefix for saved result files.
#' @param ch Logical. If \code{TRUE}, models "Change from Baseline" (\code{marker - marker_bl}).
#' @param bc Numeric (optional). Confidence level for AUC. If \code{NULL}, uses 
#'   Bonferroni correction: $1 - 0.05 / (\text{total combinations} \times 6)$.
#'
#' @return A data frame (\code{OUT}) containing:
#' \itemize{
#'   \item \code{AUC, LB, UB}: DeLong AUC and Bonferroni-corrected intervals.
#'   \item \code{Sensitivity, Specificity, PPV, NPV}: Metrics at the optimal threshold.
#'   \item \code{W_LB, W_UB}: Wilson Score Confidence Intervals for Sensitivity.
#'   \item \code{W_x, W_n}: The count of True Positives and Total Positives used for Wilson calculation.
#'   \item \code{Threshold_label/Value}: The optimal cutoff as a Probability or raw Marker value.
#' }
#'
#' @details
#' \strong{1. Combinatorial Logic:}
#' The function generates the power set of all biomarkers in \code{mar_nam}. For $n$ 
#' biomarkers, it tests $2^n - 1$ combinations per dataset.
#'
#' \strong{2. Data Processing & Modeling:}
#' For each combination, it merges endpoint data with baseline characteristics and 
#' fits **six models**: Unadjusted and Adjusted versions of Partly Conditional (PC), 
#' Time-Dependent (TD) Cox, and Simple Cox PH.
#'
#' 
#'
#' \strong{3. Performance Metrics:}
#' \itemize{
#'   \item \bold{AUC:} Calculated via \code{pROC::roc}. CIs are derived using the 
#'     DeLong method.
#'   \item \bold{Wilson Intervals:} Unlike standard Wald intervals, the Wilson 
#'     Score interval (\code{binom.confint}) is used for Sensitivity to provide 
#'     more robust coverage, especially near 0 or 1.
#'   \item \bold{Optimal Threshold:} Identified via Youden’s J. The \code{Threshold_label} 
#'     logic determines if the risk is monotonic with the biomarker to report the 
#'     cutoff in original units (e.g., pg/mL) rather than just probability.
#' }
#'
#' @import cadypredictions
#' @import tidyverse
#' @import partlyconditional
#' @import survival
#' @import pROC
#' @import binom
#' @seealso
#' The model fitting uses \code{\link[cadypredictions]{predict_risk}}.
#' The AUC and classification metrics rely heavily on the \code{pROC} package.
#'
#' @examples
#' \dontrun{
#' # Analyze BNP and CRP for a subset of endpoints, adjusting for only age and LVEF
#' custom_predictors <- c("Age", "lvef_mp_bas")
#' RESLIST_AUC(
#'   dir_in = "C:/data_location/",
#'   dir_out = "C:/results_location/",
#'   dat_nam = c("cady_data_ct", "cady_data_drug"),
#'   mar_nam = c("BNP", "CRP"),
#'   predictores = custom_predictors
#' )
#' }
#'
#' @export
#' 
#' 

RESLIST_AUC <- function(dir_in = "M:/CRF/ICORG/Studies/CADY/Clinical_Study_Report/Report/data/",
                        dir_out = "M:/CRF/ICORG/Studies/CADY/Clinical_Study_Report/Report/results/",
                        dat_nam = c("cady_data_ct","cady_data_drug","cady_data_max_50","cady_data_max_53","cady_data_mp_50","cady_data_mp_53"),
                        mar_nam = c("BNP","NT_pro_BNP","CRP","hsTnI_STAT","Galectin_3"),
                        predictores = c("Age","lvef_mp_bas","diabetes_mellitus_YN","hypertension_YN","dyslipidemia_YN","treatment_reg"),
                        saveopt = NULL,
                        ch = FALSE,
                        bc = NULL){
  
  
  library(cadypredictions)
  library(tidyverse)
  library(partlyconditional)
  library(survival)
  library(pROC)
  library(binom)
  
  # Loop from length 1 to length 4
  all_combs <- lapply(1:length(mar_nam), function(x) combn(mar_nam, x, simplify = FALSE))
  # Flatten the list
  all_combs <- unlist(all_combs, recursive = FALSE)
  numdim <- length(all_combs) * length(dat_nam)
  if (is.null(bc)) bc <- 1-0.05/(numdim*4)
  
  
  datasets_out <- lapply(dat_nam,function (dn){
    
    biomarker_out <- lapply(1:length(all_combs),function(bio){
      
      marker_name <- all_combs[[bio]] 
      marker_name_lab <- marker_name
      if (ch) marker_name_lab <- paste(marker_name," (Change)",sep="")
      
      data_name <- dn
      cat("\n\n\n Row = ",bio,"; Data and Biomarker: ",data_name," and ",marker_name_lab,"\n\n")  
      
      pa <- paste(dir_in,data_name,".csv",sep="")
      masterD <- read.csv(pa)
      n0 <- "SubjectID"
      n1 <- mar_nam
      n2 <- paste(mar_nam,"_bl",sep="")
      n3 <- c("lvef_mp_bas","lvef_max_bas","time_to_event","status","time_to_sample")
      n4 <- predictores
      wher <- names(masterD) %in% c(n0,n1,n2,n3,n4)
      masterD <- masterD[,wher]
      
      
      bldata <- read.csv(paste(dir_in,"baseline_data.csv",sep=""))
      wher1 <- names(bldata) %in% predictores
      wher2 <- names(bldata) %in% names(masterD)
      wher <- wher1 & !wher2
      
      if (!all(!wher) ) {
        bldata <- bldata[,c("SubjectID",names(bldata)[wher])]
        masterD <- masterD %>% left_join(bldata,by = "SubjectID")
      }
      
      ###################
      ### Unadjusted PC
      ###################
      fit <- predict_risk(
        datos = masterD,
        marker_name = marker_name,
        log_marker = TRUE,
        pred_from = 150,
        pred_to = 240,
        change = ch
      )
      
      df <- as.data.frame(fit$pred_data)
      
      # --- 2. ROC Calculation (pROC) ---
      what <- "Risk_PC"
      ddff <- df[,c("status",what)] %>% na.omit()
      roc_obj <- pROC::roc(ddff$status, ddff[, what], algorithm = 1, quiet = TRUE)
      
      # Calculate AUC and CI
      ci_delong <- ci(roc_obj, method = "delong",conf.level=bc)
      ci_delong_nobon <- ci(roc_obj, method = "delong",conf.level=0.95)
      
      auc <- round(ci_delong[2], 3)
      lb <- round(ci_delong[1], 3)
      ub <- round(ci_delong[3], 3)
      lb_nobon <- round(ci_delong_nobon[1], 3)
      ub_nobon <- round(ci_delong_nobon[3], 3)
      
      # --- 3. Optimal Threshold Calculation (Youden's J) ---
      optimal_coords <- coords(roc_obj, x = "best", best.method = "youden",
                               ret = c("threshold", "specificity", "sensitivity", "ppv", "npv","tp","fn","fp","tn"))
      if (dim(optimal_coords)[1]>1) {
        optimal_coords <- optimal_coords %>% arrange(threshold) %>% slice(1)
      }
      
      # 2. Calculate Wilson Interval Sensitivity = TP / (TP + FN)
      wci <- binom.confint(x = optimal_coords$tp, 
                           n = (optimal_coords$tp + optimal_coords$fn), 
                           methods = "wilson")
      W_x = wci$x
      W_n = wci$n
      W_LB = wci$lower
      W_UB = wci$upper
      
      sen <- round(optimal_coords$sensitivity,3)
      spe <- round(optimal_coords$specificity,3)
      ppv <- round(optimal_coords$ppv,3)
      npv <- round(optimal_coords$npv,3)
      
      
      # Conditional Optimal Threshold Label Logic
      df$dup <- df[, what]
      temp1 <- order(df$marker1)
      temp2 <- order(df$dup)
      if (all(temp1 == temp2)) {
        # If ordering is the same, use the marker value
        df1 <- df %>% 
          arrange(dup) %>% 
          filter(dup <= optimal_coords$threshold) %>% 
          slice_tail(n = 1)
        threshold_label <- "Marker"
        threshold_value <- round(df1$marker1, 1)
      } else {
        # If ordering is different, use the raw probability threshold
        threshold_label <-"Prob" 
        threshold_value <- round(optimal_coords$threshold, 3)
      }
      
      out1 <- data.frame(data_name=data_name,
                         type = "Unadjusted",
                         model = "Partly conditional",
                         AUC = auc,
                         LB = lb,
                         UB = ub,
                         LB_nobon = lb_nobon,
                         UB_nobon = ub_nobon,
                         Threshold_label = threshold_label,
                         Threshold_value = threshold_value,
                         Sensitivity = sen,
                         Specificity = spe,
                         PPV = ppv,
                         NPV = npv,
                         W_x = W_x,
                         W_n = W_n,
                         W_LB = W_LB,
                         W_UB = W_UB,
                         PRED = paste(fit$pc_model_predictors,collapse="+"))
      
      
      ###################
      ### Unadjusted TD
      ###################
      df <- as.data.frame(fit$pred_data)
      
      # --- 2. ROC Calculation (pROC) ---
      what <- "Risk_TDcox"
      ddff <- df[,c("status",what)] %>% na.omit()
      roc_obj <- pROC::roc(ddff$status, ddff[, what], algorithm = 1, quiet = TRUE)
      
      
      
      # Calculate AUC and CI
      ci_delong <- ci(roc_obj, method = "delong",conf.level=bc)
      ci_delong_nobon <- ci(roc_obj, method = "delong",conf.level=0.95)
      
      auc <- round(ci_delong[2], 3)
      lb <- round(ci_delong[1], 3)
      ub <- round(ci_delong[3], 3)
      lb_nobon <- round(ci_delong_nobon[1], 3)
      ub_nobon <- round(ci_delong_nobon[3], 3)
      
      # --- 3. Optimal Threshold Calculation (Youden's J) ---
      optimal_coords <- coords(roc_obj, x = "best", best.method = "youden",
                               ret = c("threshold", "specificity", "sensitivity", "ppv", "npv","tp","fn","fp","tn"))
      if (dim(optimal_coords)[1]>1) {
        optimal_coords <- optimal_coords %>% arrange(threshold) %>% slice(1)
      }
      
      # 2. Calculate Wilson Interval Sensitivity = TP / (TP + FN)
      wci <- binom.confint(x = optimal_coords$tp, 
                           n = (optimal_coords$tp + optimal_coords$fn), 
                           methods = "wilson")
      W_x = wci$x
      W_n = wci$n
      W_LB = wci$lower
      W_UB = wci$upper
      
      sen <- round(optimal_coords$sensitivity,3)
      spe <- round(optimal_coords$specificity,3)
      ppv <- round(optimal_coords$ppv,3)
      npv <- round(optimal_coords$npv,3)
      
      
      # Conditional Optimal Threshold Label Logic
      df$dup <- df[, what]
      temp1 <- order(df$marker1)
      temp2 <- order(df$dup)
      if (all(temp1 == temp2)) {
        # If ordering is the same, use the marker value
        df1 <- df %>% 
          arrange(dup) %>% 
          filter(dup <= optimal_coords$threshold) %>% 
          slice_tail(n = 1)
        threshold_label <- "Marker"
        threshold_value <- round(df1$marker1, 1)
      } else {
        # If ordering is different, use the raw probability threshold
        threshold_label <-"Prob" 
        threshold_value <- round(optimal_coords$threshold, 3)
      }
      
      out2 <- data.frame(data_name=data_name,
                         type = "Unadjusted",
                         model = "Time dependent Cox PH",
                         AUC = auc,
                         LB = lb,
                         UB = ub,
                         LB_nobon = lb_nobon,
                         UB_nobon = ub_nobon,
                         Threshold_label = threshold_label,
                         Threshold_value = threshold_value,
                         Sensitivity = sen,
                         Specificity = spe,
                         PPV = ppv,
                         NPV = npv,
                         W_x = W_x,
                         W_n = W_n,
                         W_LB = W_LB,
                         W_UB = W_UB,
                         PRED = paste(fit$tdcox_model_predictors,collapse="+"))
      
      
      ##########################
      ### Unadjusted Simple Cox
      ##########################
      df <- as.data.frame(fit$pred_data)
      
      # --- 2. ROC Calculation (pROC) ---
      what <- "Risk_cox_simple"
      ddff <- df[,c("status",what)] %>% na.omit()
      roc_obj <- pROC::roc(ddff$status, ddff[, what], algorithm = 1, quiet = TRUE)
      
      # Calculate AUC and CI
      ci_delong <- ci(roc_obj, method = "delong",conf.level=bc)
      ci_delong_nobon <- ci(roc_obj, method = "delong",conf.level=0.95)
      
      auc <- round(ci_delong[2], 3)
      lb <- round(ci_delong[1], 3)
      ub <- round(ci_delong[3], 3)
      lb_nobon <- round(ci_delong_nobon[1], 3)
      ub_nobon <- round(ci_delong_nobon[3], 3)
      
      # --- 3. Optimal Threshold Calculation (Youden's J) ---
      optimal_coords <- coords(roc_obj, x = "best", best.method = "youden",
                               ret = c("threshold", "specificity", "sensitivity", "ppv", "npv","tp","fn","fp","tn"))
      if (dim(optimal_coords)[1]>1) {
        optimal_coords <- optimal_coords %>% arrange(threshold) %>% slice(1)
      }
      
      # 2. Calculate Wilson Interval Sensitivity = TP / (TP + FN)
      wci <- binom.confint(x = optimal_coords$tp, 
                           n = (optimal_coords$tp + optimal_coords$fn), 
                           methods = "wilson")
      W_x = wci$x
      W_n = wci$n
      W_LB = wci$lower
      W_UB = wci$upper
      
      
      sen <- round(optimal_coords$sensitivity,3)
      spe <- round(optimal_coords$specificity,3)
      ppv <- round(optimal_coords$ppv,3)
      npv <- round(optimal_coords$npv,3)
      
      
      # Conditional Optimal Threshold Label Logic
      df$dup <- df[, what]
      temp1 <- order(df$marker1)
      temp2 <- order(df$dup)
      if (all(temp1 == temp2)) {
        # If ordering is the same, use the marker value
        df1 <- df %>% 
          arrange(dup) %>% 
          filter(dup <= optimal_coords$threshold) %>% 
          slice_tail(n = 1)
        threshold_label <- "Marker"
        threshold_value <- round(df1$marker1, 1)
      } else {
        # If ordering is different, use the raw probability threshold
        threshold_label <-"Prob" 
        threshold_value <- round(optimal_coords$threshold, 3)
      }
      
      out3 <- data.frame(data_name=data_name,
                         type = "Unadjusted",
                         model = "Simple Cox PH",
                         AUC = auc,
                         LB = lb,
                         UB = ub,
                         LB_nobon = lb_nobon,
                         UB_nobon = ub_nobon,
                         Threshold_label = threshold_label,
                         Threshold_value = threshold_value,
                         Sensitivity = sen,
                         Specificity = spe,
                         PPV = ppv,
                         NPV = npv,
                         W_x = W_x,
                         W_n = W_n,
                         W_LB = W_LB,
                         W_UB = W_UB,
                         PRED = paste(fit$simplecox_model_predictors,collapse="+"))
      
      
      ################
      ### Adjusted PC
      ################
      fit <- predict_risk(
        datos = masterD,
        marker_name = marker_name,
        Predictors = predictores,
        log_marker = TRUE,
        pred_from = 150,
        pred_to = 240,
        change = ch
      )
      pc_model_predictors<-fit$pc_model_predictors
      tdcox_model_predictors <- fit$tdcox_model_predictors
      simplecox_model_predictors <- fit$simplecox_model_predictors
      
      
      df <- as.data.frame(fit$pred_data)
      
      # --- 2. ROC Calculation (pROC) ---
      what = "Risk_PC"
      ddff <- df[,c("status",what)] %>% na.omit()
      roc_obj <- pROC::roc(ddff$status, ddff[, what], algorithm = 1, quiet = TRUE)
      
      # Calculate AUC and CI
      ci_delong <- ci(roc_obj, method = "delong",conf.level=bc)
      ci_delong_nobon <- ci(roc_obj, method = "delong",conf.level=0.95)
      
      auc <- round(ci_delong[2], 3)
      lb <- round(ci_delong[1], 3)
      ub <- round(ci_delong[3], 3)
      lb_nobon <- round(ci_delong_nobon[1], 3)
      ub_nobon <- round(ci_delong_nobon[3], 3)
      
      # --- 3. Optimal Threshold Calculation (Youden's J) ---
      optimal_coords <- coords(roc_obj, x = "best", best.method = "youden",
                               ret = c("threshold", "specificity", "sensitivity", "ppv", "npv","tp","fn","fp","tn"))
      if (dim(optimal_coords)[1]>1) {
        optimal_coords <- optimal_coords %>% arrange(threshold) %>% slice(1)
      }
      
      
      # 2. Calculate Wilson Interval Sensitivity = TP / (TP + FN)
      wci <- binom.confint(x = optimal_coords$tp, 
                           n = (optimal_coords$tp + optimal_coords$fn), 
                           methods = "wilson")
      W_x = wci$x
      W_n = wci$n
      W_LB = wci$lower
      W_UB = wci$upper
      
      
      sen <- round(optimal_coords$sensitivity,3)
      spe <- round(optimal_coords$specificity,3)
      ppv <- round(optimal_coords$ppv,3)
      npv <- round(optimal_coords$npv,3)
      
      
      # Conditional Optimal Threshold Label Logic
      df$dup <- df[, what]
      temp1 <- order(df$marker1)
      temp2 <- order(df$dup)
      if (all(temp1 == temp2)) {
        # If ordering is the same, use the marker value
        df1 <- df %>% 
          arrange(dup) %>% 
          filter(dup <= optimal_coords$threshold) %>% 
          slice_tail(n = 1)
        threshold_label <- "Marker"
        threshold_value <- round(df1$marker1, 1)
      } else {
        # If ordering is different, use the raw probability threshold
        threshold_label <-"Prob" 
        threshold_value <- round(optimal_coords$threshold, 3)
      }
      
      out4 <- data.frame(data_name=data_name,
                         type = "Adjusted",
                         model = "Partly conditional",
                         AUC = auc,
                         LB = lb,
                         UB = ub,
                         LB_nobon = lb_nobon,
                         UB_nobon = ub_nobon,
                         Threshold_label = threshold_label,
                         Threshold_value = threshold_value,
                         Sensitivity = sen,
                         Specificity = spe,
                         PPV = ppv,
                         NPV = npv,
                         W_x = W_x,
                         W_n = W_n,
                         W_LB = W_LB,
                         W_UB = W_UB,
                         PRED = paste(fit$pc_model_predictors,collapse="+"))
      
      
      ###################
      ### Adjusted TD
      ###################
      # --- 2. ROC Calculation (pROC) ---
      what = "Risk_TDcox"
      ddff <- df[,c("status",what)] %>% na.omit()
      roc_obj <- pROC::roc(ddff$status, ddff[, what], algorithm = 1, quiet = TRUE)
      
      
      # Calculate AUC and CI
      ci_delong <- ci(roc_obj, method = "delong",conf.level=bc)
      ci_delong_nobon <- ci(roc_obj, method = "delong",conf.level=0.95)
      
      auc <- round(ci_delong[2], 3)
      lb <- round(ci_delong[1], 3)
      ub <- round(ci_delong[3], 3)
      lb_nobon <- round(ci_delong_nobon[1], 3)
      ub_nobon <- round(ci_delong_nobon[3], 3)
      
      # --- 3. Optimal Threshold Calculation (Youden's J) ---
      optimal_coords <- coords(roc_obj, x = "best", best.method = "youden",
                               ret = c("threshold", "specificity", "sensitivity", "ppv", "npv","tp","fn","fp","tn"))
      if (dim(optimal_coords)[1]>1) {
        optimal_coords <- optimal_coords %>% arrange(threshold) %>% slice(1)
      }
      
      # 2. Calculate Wilson Interval Sensitivity = TP / (TP + FN)
      wci <- binom.confint(x = optimal_coords$tp, 
                           n = (optimal_coords$tp + optimal_coords$fn), 
                           methods = "wilson")
      W_x = wci$x
      W_n = wci$n
      W_LB = wci$lower
      W_UB = wci$upper
      
      
      sen <- round(optimal_coords$sensitivity,3)
      spe <- round(optimal_coords$specificity,3)
      ppv <- round(optimal_coords$ppv,3)
      npv <- round(optimal_coords$npv,3)
      
      
      # Conditional Optimal Threshold Label Logic
      df$dup <- df[, what]
      temp1 <- order(df$marker1)
      temp2 <- order(df$dup)
      if (all(temp1 == temp2)) {
        # If ordering is the same, use the marker value
        df1 <- df %>% 
          arrange(dup) %>% 
          filter(dup <= optimal_coords$threshold) %>% 
          slice_tail(n = 1)
        threshold_label <- "Marker"
        threshold_value <- round(df1$marker1, 1)
      } else {
        # If ordering is different, use the raw probability threshold
        threshold_label <-"Prob" 
        threshold_value <- round(optimal_coords$threshold, 3)
      }
      
      out5 <- data.frame(data_name=data_name,
                         type = "Adjusted",
                         model = "Time dependent Cox PH",
                         AUC = auc,
                         LB = lb,
                         UB = ub,
                         LB_nobon = lb_nobon,
                         UB_nobon = ub_nobon,
                         Threshold_label = threshold_label,
                         Threshold_value = threshold_value,
                         Sensitivity = sen,
                         Specificity = spe,
                         PPV = ppv,
                         NPV = npv,
                         W_x = W_x,
                         W_n = W_n,
                         W_LB = W_LB,
                         W_UB = W_UB,
                         PRED = paste(fit$tdcox_model_predictors,collapse="+"))
      
      
      ########################
      ### Adjusted Simple Cox
      ########################
      # --- 2. ROC Calculation (pROC) ---
      what = "Risk_cox_simple"
      ddff <- df[,c("status",what)] %>% na.omit()
      roc_obj <- pROC::roc(ddff$status, ddff[, what], algorithm = 1, quiet = TRUE)
      
      # Calculate AUC and CI
      ci_delong <- ci(roc_obj, method = "delong",conf.level=bc)
      ci_delong_nobon <- ci(roc_obj, method = "delong",conf.level=0.95)
      
      auc <- round(ci_delong[2], 3)
      lb <- round(ci_delong[1], 3)
      ub <- round(ci_delong[3], 3)
      lb_nobon <- round(ci_delong_nobon[1], 3)
      ub_nobon <- round(ci_delong_nobon[3], 3)
      
      # --- 3. Optimal Threshold Calculation (Youden's J) ---
      optimal_coords <- coords(roc_obj, x = "best", best.method = "youden",
                               ret = c("threshold", "specificity", "sensitivity", "ppv", "npv","tp","fn","fp","tn"))
      if (dim(optimal_coords)[1]>1) {
        optimal_coords <- optimal_coords %>% arrange(threshold) %>% slice(1)
      }
      
      # 2. Calculate Wilson Interval Sensitivity = TP / (TP + FN)
      wci <- binom.confint(x = optimal_coords$tp, 
                                 n = (optimal_coords$tp + optimal_coords$fn), 
                                 methods = "wilson")
      W_x = wci$x
      W_n = wci$n
      W_LB = wci$lower
      W_UB = wci$upper
      
      sen <- round(optimal_coords$sensitivity,3)
      spe <- round(optimal_coords$specificity,3)
      ppv <- round(optimal_coords$ppv,3)
      npv <- round(optimal_coords$npv,3)
      
      
      # Conditional Optimal Threshold Label Logic
      df$dup <- df[, what]
      temp1 <- order(df$marker1)
      temp2 <- order(df$dup)
      if (all(temp1 == temp2)) {
        # If ordering is the same, use the marker value
        df1 <- df %>% 
          arrange(dup) %>% 
          filter(dup <= optimal_coords$threshold) %>% 
          slice_tail(n = 1)
        threshold_label <- "Marker"
        threshold_value <- round(df1$marker1, 1)
      } else {
        # If ordering is different, use the raw probability threshold
        threshold_label <-"Prob" 
        threshold_value <- round(optimal_coords$threshold, 3)
      }
      
      out6 <- data.frame(data_name=data_name,
                         type = "Adjusted",
                         model = "Simple Cox PH",
                         AUC = auc,
                         LB = lb,
                         UB = ub,
                         LB_nobon = lb_nobon,
                         UB_nobon = ub_nobon,
                         Threshold_label = threshold_label,
                         Threshold_value = threshold_value,
                         Sensitivity = sen,
                         Specificity = spe,
                         PPV = ppv,
                         NPV = npv,
                         W_x = W_x,
                         W_n = W_n,
                         W_LB = W_LB,
                         W_UB = W_UB,
                         PRED = paste(fit$simplecox_model_predictors,collapse="+"))
      
      
      
      
      out <- rbind(out1,out2,out3,out4,out5,out6)
      
      row.names(out) <- NULL
      mb <- paste(marker_name_lab,collapse=" + ")
      out <- out %>% mutate(model_bio = mb)
      out <- out %>% select(data_name, model_bio,type:PRED)
      
      out
      
    })
    
    biomarker_out <- do.call("rbind",biomarker_out)
    
    
    biomarker_out
  })
  
  
  OUT <- do.call("rbind",datasets_out)
  
  if (!is.null(saveopt)) {
    save(OUT,file = paste(dir_out,saveopt,".RData",sep="") )
    write.csv(OUT,file = paste(dir_out,saveopt,".csv",sep=""),row.names = FALSE)
  }
  
  OUT
  
  
}



