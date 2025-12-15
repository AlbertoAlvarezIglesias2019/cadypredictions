#' Run Survival Models Across Multiple Biomarkers and Clinical Endpoints
#'
#' This function systematically loads various datasets (corresponding to different
#' clinical endpoints or definitions, e.g., CT definitions), iterates through a
#' set of predefined biomarkers, and fits both unadjusted and adjusted survival
#' models. The function utilizes both the Partly Conditional (PC) model and
#' the standard Time-Dependent (TD) Cox Proportional Hazards model for risk
#' prediction.
#'
#' @param dir_in A character string specifying the path to the directory
#'   where the input CSV data files are located. Defaults to a specific network
#'   path. The function expects files named like \code{cady_data_ct.csv} and
#'   \code{baseline_data.csv}.
#' @param dir_out A character string specifying the path to the directory
#'   where the final results (\code{RESULTS.RData} and \code{RESULTS.csv})
#'   should be saved. Defaults to a specific network path.
#'
#' @return The function does not explicitly return a value but saves the
#'   aggregated results as an RData file and a CSV file named \code{RESULTS.RData}
#'   and \code{RESULTS.csv} in the directory specified by \code{dir_out}.
#'   The saved data frame contains the Hazard Ratios (HR), Confidence Intervals
#'   (LB/UB), and P-values (PVAL) for every combination of marker, data definition,
#'   adjustment status, and model type.
#'
#' @details
#' The function performs a comprehensive survival analysis for risk prediction:
#'
#' \enumerate{
#'   \item **Setup:** Defines a grid (\code{pointer}) covering all $6 \times 5 = 30$
#'     combinations of data definitions (e.g., CT outcome types) and biomarkers.
#'   \item **Data Loading & Merging:** For each combination:
#'     \itemize{
#'       \item Loads the main outcome data (e.g., \code{cady_data_ct.csv}).
#'       \item Loads baseline covariate data (\code{baseline_data.csv}).
#'       \item Merges the two data frames by \code{SubjectID}.
#'     \itemize}
#'   \item **Model Fitting (using \code{predict_risk} from \code{cadypredictions}):**
#'     Four models are fitted for each combination of data and marker:
#'     \itemize{
#'       \item **Unadjusted Partly Conditional (PC):** Includes only the log-transformed
#'         biomarker.
#'       \item **Unadjusted Time-Dependent (TD) Cox PH:** Includes only the log-transformed
#'         biomarker.
#'       \item **Adjusted Partly Conditional (PC):** Includes the log-transformed biomarker
#'         and a set of baseline covariates (\code{Age}, \code{lvef_mp_bas},
#'         \code{diabetes_mellitus_YN}, \code{hypertension_YN}, \code{dyslipidemia_YN},
#'         \code{treatment_reg}).
#'       \item **Adjusted Time-Dependent (TD) Cox PH:** Includes the log-transformed
#'         biomarker and the same baseline covariates.
#'     \itemize}
#'   \item **Result Extraction:** Hazard Ratios (HR), 95\% Confidence Intervals (LB, UB),
#'     and p-values (PVAL) are extracted for the biomarker term (\code{marker1})
#'     from the model output and rounded.
#'   \item **Aggregation:** The results from all $4 \times 30 = 120$ models are
#'     combined into a single data frame (\code{OUT}).
#'   \item **Saving:** The final results table is saved as both an RData object
#'     and a CSV file.
#' }
#'
#' @import cadypredictions
#' @import tidyverse
#' @import partlyconditional
#' @import survival
#' @importFrom base paste read.csv data.frame round format
#' @importFrom dplyr select left_join
#' @importFrom stats confint
#'
#' @seealso
#' The core modeling is performed by \code{\link[cadypredictions]{predict_risk}}.
#'
#' @examples
#' \dontrun{
#' # Assuming the necessary data files are in a local temporary directory:
#' # e.g., 'C:/my_project/data/' and results go to 'C:/my_project/results/'
#' # WARNING: This function depends on the existence of specific files
#' # and the 'cadypredictions' and 'partlyconditional' packages.
#'
#' RESULTS(
#'   dir_in = "C:/my_project/data/",
#'   dir_out = "C:/my_project/results/"
#' )
#' # The results will be saved in "C:/my_project/results/RESULTS.csv"
#' }
#'
#' @export
#' 
#' 

RESULTS <- function(dir_in = "M:/CRF/ICORG/Studies/CADY/Clinical_Study_Report/Report/data/",
                   dir_out = "M:/CRF/ICORG/Studies/CADY/Clinical_Study_Report/Report/results/"){
  
  library(cadypredictions)
  library(tidyverse)
  library(partlyconditional)
  library(survival)
  
  
  pointer <- expand.grid(data_name = c("cady_data_ct","cady_data_drug","cady_data_max_50","cady_data_max_53","cady_data_mp_50","cady_data_mp_53"),
                         marker_name = c("BNP","NT_pro_BNP","CRP","hsTnI_STAT","Galectin_3"))
  
  temp <- lapply(1:dim(pointer)[1],function (kkk){
    #data_name <- "cady_data_mp_50"
    #marker_name = "NT_pro_BNP"
    data_name <- as.character(pointer$data_name[kkk])
    marker_name <- as.character(pointer$marker_name[kkk])
    
    #data_name <- "cady_data_ct" 
    #marker_name <- "hsTnI_STAT"
    
    cat("\n\n\n",data_name," and ",marker_name,"\n\n")  
    pa <- paste(dir_in,data_name,".csv",sep="")
    masterD <- read.csv(pa)
    masterD <- masterD %>% select(SubjectID,BNP,NT_pro_BNP,CRP,hsTnI_STAT,Galectin_3,
                                  BNP_bl,NT_pro_BNP_bl,CRP_bl,hsTnI_STAT_bl,Galectin_3_bl,
                                  lvef_mp_bas,lvef_max_bas,time_to_event,status,time_to_sample)
    
    bldata <- read.csv(paste(dir_in,"baseline_data.csv",sep=""))
    bldata <- bldata %>% select(SubjectID,Age,diabetes_mellitus_YN,hypertension_YN,
                                dyslipidemia_YN,treatment_reg)
    
    masterD <- masterD %>% left_join(bldata,by = "SubjectID")
    
    
    ###################
    ### Unadjusted PC
    ###################
    fit <- predict_risk(
      datos = masterD,
      marker_name = marker_name,
      log_marker = TRUE,
      pred_from = 150,
      pred_to = 240
    )
    
    hr <- exp(as.numeric(fit$pc_model$coefficients["marker1"]))
    wher <- row.names(confint(fit$pc_model)) %in% "marker1"
    lb <- exp(as.numeric(confint(fit$pc_model)[wher,][1]))
    ub <- exp(as.numeric(confint(fit$pc_model)[wher,][2]))
    hr <- round(hr,2)
    lb <- round(lb,2)
    ub <- round(ub,2)
    wher <- row.names(summary(fit$pc_model)$coefficients) %in% "marker1"
    pval <- summary(fit$pc_model)$coefficients[wher,"Pr(>|z|)"]
    pval <- round(pval,3)
    pval <- format(pval,nsmall=3)
    
    out1 <- data.frame(marker_name=marker_name,
                       data_name=data_name,
                       type = "Unadjusted",
                       model = "Partly conditional",
                       HR = hr,
                       LB = lb,
                       UB = ub,
                       PVAL = pval)
    
    ###################
    ### Unadjusted TD
    ###################
    hr <- exp(as.numeric(fit$tdcox_model$coefficients["marker1"]))
    wher <- row.names(confint(fit$tdcox_model)) %in% "marker1"
    lb <- exp(as.numeric(confint(fit$tdcox_model)[wher,][1]))
    ub <- exp(as.numeric(confint(fit$tdcox_model)[wher,][2]))
    hr <- round(hr,2)
    lb <- round(lb,2)
    ub <- round(ub,2)
    wher <- row.names(summary(fit$tdcox_model)$coefficients) %in% "marker1"
    pval <- summary(fit$tdcox_model)$coefficients[wher,"Pr(>|z|)"]
    pval <- round(pval,3)
    pval <- format(pval,nsmall=3)
    
    out2 <- data.frame(marker_name=marker_name,
                       data_name=data_name,
                       type = "Unadjusted",
                       model = "Time dependent Cox PH",
                       HR = hr,
                       LB = lb,
                       UB = ub,
                       PVAL = pval)
    
    
    ################
    ### Adjusted PC
    ################
    fit <- predict_risk(
      datos = masterD,
      marker_name = marker_name,
      Predictors = c("Age","lvef_mp_bas","diabetes_mellitus_YN","hypertension_YN","dyslipidemia_YN","treatment_reg"),
      log_marker = TRUE,
      pred_from = 150,
      pred_to = 240
    )
    
    hr <- exp(as.numeric(fit$pc_model$coefficients["marker1"]))
    wher <- row.names(confint(fit$pc_model)) %in% "marker1"
    lb <- exp(as.numeric(confint(fit$pc_model)[wher,][1]))
    ub <- exp(as.numeric(confint(fit$pc_model)[wher,][2]))
    hr <- round(hr,2)
    lb <- round(lb,2)
    ub <- round(ub,2)
    wher <- row.names(summary(fit$pc_model)$coefficients) %in% "marker1"
    pval <- summary(fit$pc_model)$coefficients[wher,"Pr(>|z|)"]
    pval <- round(pval,3)
    pval <- format(pval,nsmall=3)
    
    out3 <- data.frame(marker_name=marker_name,
                       data_name=data_name,
                       type = "Adjusted",
                       model = "Partly conditional",
                       HR = hr,
                       LB = lb,
                       UB = ub,
                       PVAL = pval)
    
    
    
    
    
    ################
    ### Adjusted TD
    ################
    hr <- exp(as.numeric(fit$tdcox_model$coefficients["marker1"]))
    wher <- row.names(confint(fit$tdcox_model)) %in% "marker1"
    lb <- exp(as.numeric(confint(fit$tdcox_model)[wher,][1]))
    ub <- exp(as.numeric(confint(fit$tdcox_model)[wher,][2]))
    hr <- round(hr,2)
    lb <- round(lb,2)
    ub <- round(ub,2)
    wher <- row.names(summary(fit$tdcox_model)$coefficients) %in% "marker1"
    pval <- summary(fit$tdcox_model)$coefficients[wher,"Pr(>|z|)"]
    pval <- round(pval,3)
    pval <- format(pval,nsmall=3)
    
    out4 <- data.frame(marker_name=marker_name,
                       data_name=data_name,
                       type = "Adjusted",
                       model = "Time dependent Cox PH",
                       HR = hr,
                       LB = lb,
                       UB = ub,
                       PVAL = pval)
    
    
    out <- rbind(out1,out2,out3,out4)
    out
  })
  
  OUT <- do.call("rbind",temp)
  
  save(OUT,file = paste(dir_out,"RESULTS.RData",sep="") )
  write.csv(OUT,file = paste(dir_out,"RESULTS.csv",sep=""),row.names = FALSE)
}


