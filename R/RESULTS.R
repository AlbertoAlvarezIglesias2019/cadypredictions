#' Run Survival Models Across Multiple Biomarkers and Clinical Endpoints
#'
#' This function systematically loads various datasets (corresponding to different
#' clinical endpoints), iterates through a user-defined set of biomarkers, and fits
#' both unadjusted and *adjusted* survival models. The analysis uses both the
#' Partly Conditional (PC) model and the standard Time-Dependent (TD) Cox
#' Proportional Hazards model to estimate Hazard Ratios (HRs).
#'
#' @param dir_in A character string specifying the path to the directory
#'   where the input CSV data files are located. Defaults to a specific network
#'   path. The function expects files named like \code{dat_nam[i]}.csv and
#'   \code{baseline_data.csv}.
#' @param dir_out A character string specifying the path to the directory
#'   where the final results (\code{RESULTS.RData} and \code{RESULTS.csv})
#'   should be saved. Defaults to a specific network path.
#' @param dat_nam A character vector listing the base names of the data files
#'   representing different clinical endpoint definitions (e.g., CT outcome types).
#'   Defaults to a predefined set of six data names.
#' @param mar_nam A character vector listing the biomarker names to analyze.
#'   Must correspond to columns in the loaded data. Defaults to a predefined set of five biomarkers.
#' @param predictores A character vector specifying the names of the covariates
#'   to be included in the 'Adjusted' models. These variables must be present
#'   in the \code{baseline_data.csv} file. Defaults to a predefined set of six covariates.
#' @param saveopt A character string (optional). If provided, the final aggregated
#'   results data frame (\code{OUT}) will be saved to \code{dir_out} as a CSV and RData
#'   file using this string as the file name prefix. If \code{NULL} (default), no file is saved.
#'   
#'
#' @return A data frame (\code{OUT}) containing the summarized results for every
#'   combination of data, marker, adjustment type (Unadjusted/Adjusted), and model,
#'   including the following columns:
#' \itemize{
#'   \item \code{marker_name}: The biomarker analyzed.
#'   \item \code{data_name}: The dataset used.
#'   \item \code{type}: "Unadjusted" or "Adjusted".
#'   \item \code{model}: The type of survival model ("Partly conditional", "Time dependent Cox PH", or "Simple Cox PH").
#'   \item \code{HR}: The calculated Hazard Ratio for the biomarker ($\text{log}_2$ transformed). 
#'   \item \code{LB}, \code{UB}: The lower and upper bounds of the 95% Confidence Interval for the HR.
#'   \item \code{PVAL}: The p-value for the biomarker's coefficient.
#'   \item \code{PRED}: A string listing the final set of predictors used in the model (this will dynamically change if the model fitting required dropping predictors due to convergence warnings).
#' }
#'
#' @details
#' The function first defines a grid (\code{pointer}) encompassing all combinations
#' of the provided \code{dat_nam} and \code{mar_nam}. It then iterates through this grid:
#'
#' \enumerate{
#'   \item **Data Processing:**
#'     \itemize{
#'       \item Loads the specific endpoint data (e.g., \code{cady\_data\_ct.csv}) and baseline data (\code{baseline\_data.csv}).
#'       \item Merges data frames by \code{SubjectID}.
#'     \itemize}
#'   \item **Model Fitting (using \code{predict_risk}):** Four survival models are fitted for each combination of data and marker:
#'     \itemize{
#'       \item **Unadjusted Partly Conditional (PC)**: Includes only the log-transformed biomarker.
#'       \item **Unadjusted Time-Dependent (TD) Cox PH** : Includes only the log-transformed biomarker.
#'       \item **Adjusted Partly Conditional (PC)**: Includes the log-transformed biomarker and the covariates specified in the **`predictores`** argument.
#'       \item **Adjusted Time-Dependent (TD) Cox PH**: Includes the log-transformed biomarker and the covariates specified in the **`predictores`** argument.
#'     \itemize}
#'   \item **Result Extraction:** The Hazard Ratio ($\text{HR} = e^{\beta}$), its 95\% Confidence Intervals (LB, UB), and the p-value (PVAL) are extracted for the biomarker term (\code{marker1}) from the model summaries. The $\text{HR}$ and CIs are rounded to two decimal places, and the P-value to three.
#'   \item **Aggregation and Saving:** Results from all models are compiled into a single data frame (\code{OUT}) and saved to the output directory.
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
#' # Run the analysis for two specific data sets, one marker, and custom predictors
#' custom_predictors <- c("Age", "treatment_reg")
#' RESULTS(
#'   dir_in = "C:/my_data/",
#'   dir_out = "C:/my_results/",
#'   dat_nam = c("cady_data_mp_50", "cady_data_mp_53"),
#'   mar_nam = "NT_pro_BNP",
#'   predictores = custom_predictors
#' )
#' # The results will be saved in "C:/my_results/RESULTS.csv"
#' }
#'
#' @export
#' 

RESULTS <- function(dir_in = "M:/CRF/ICORG/Studies/CADY/Clinical_Study_Report/Report/data/",
                   dir_out = "M:/CRF/ICORG/Studies/CADY/Clinical_Study_Report/Report/results/",
                   dat_nam = c("cady_data_ct","cady_data_drug","cady_data_max_50","cady_data_max_53","cady_data_mp_50","cady_data_mp_53"),
                   mar_nam = c("BNP","NT_pro_BNP","CRP","hsTnI_STAT","Galectin_3"),
                   predictores = c("Age","lvef_mp_bas","diabetes_mellitus_YN","hypertension_YN","dyslipidemia_YN","treatment_reg"),
                   saveopt = NULL){
  
  library(cadypredictions)
  library(tidyverse)
  library(partlyconditional)
  library(survival)
  
  
  pointer <- expand.grid(data_name = dat_nam,
                         marker_name = mar_nam)
  
  temp <- lapply(1:dim(pointer)[1],function (kkk){
    #data_name <- "cady_data_mp_50"
    #marker_name = "NT_pro_BNP"
    data_name <- as.character(pointer$data_name[kkk])
    marker_name <- as.character(pointer$marker_name[kkk])
    
    #data_name <- "cady_data_ct" 
    #marker_name <- "hsTnI_STAT"
    
    cat("\n\n\n Row = ",kkk,"; Data and Biomarker: ",data_name," and ",marker_name,"\n\n")  
    pa <- paste(dir_in,data_name,".csv",sep="")
    masterD <- read.csv(pa)
    #masterD <- masterD %>% select(SubjectID,BNP,NT_pro_BNP,CRP,hsTnI_STAT,Galectin_3,
    #                              BNP_bl,NT_pro_BNP_bl,CRP_bl,hsTnI_STAT_bl,Galectin_3_bl,
    #                              lvef_mp_bas,lvef_max_bas,time_to_event,status,time_to_sample)
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
    pred <- paste(fit$pc_model_predictors,collapse="+")
    
    out1 <- data.frame(marker_name=marker_name,
                       data_name=data_name,
                       type = "Unadjusted",
                       model = "Partly conditional",
                       HR = hr,
                       LB = lb,
                       UB = ub,
                       PVAL = pval,
                       PRED = pred)
    
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
    pred <- paste(fit$tdcox_model_predictors,collapse="+")
    
    out2 <- data.frame(marker_name=marker_name,
                       data_name=data_name,
                       type = "Unadjusted",
                       model = "Time dependent Cox PH",
                       HR = hr,
                       LB = lb,
                       UB = ub,
                       PVAL = pval,
                       PRED = pred)
    
    
    ##########################
    ### Unadjusted Simple Cox
    ##########################
    hr <- exp(as.numeric(fit$simplecox_model$coefficients["marker1"]))
    wher <- row.names(confint(fit$simplecox_model)) %in% "marker1"
    lb <- exp(as.numeric(confint(fit$simplecox_model)[wher,][1]))
    ub <- exp(as.numeric(confint(fit$simplecox_model)[wher,][2]))
    hr <- round(hr,2)
    lb <- round(lb,2)
    ub <- round(ub,2)
    wher <- row.names(summary(fit$simplecox_model)$coefficients) %in% "marker1"
    pval <- summary(fit$simplecox_model)$coefficients[wher,"Pr(>|z|)"]
    pval <- round(pval,3)
    pval <- format(pval,nsmall=3)
    pred <- paste(fit$simplecox_model_predictors,collapse="+")
    
    out3 <- data.frame(marker_name=marker_name,
                       data_name=data_name,
                       type = "Unadjusted",
                       model = "Simple Cox PH",
                       HR = hr,
                       LB = lb,
                       UB = ub,
                       PVAL = pval,
                       PRED = pred)
    
    
    ################
    ### Adjusted PC
    ################
    fit <- predict_risk(
      datos = masterD,
      marker_name = marker_name,
      Predictors = predictores,
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
    pred <- paste(fit$pc_model_predictors,collapse="+")
    
    out4 <- data.frame(marker_name=marker_name,
                       data_name=data_name,
                       type = "Adjusted",
                       model = "Partly conditional",
                       HR = hr,
                       LB = lb,
                       UB = ub,
                       PVAL = pval,
                       PRED = pred)
    
    
    
    
    
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
    pred <- paste(fit$tdcox_model_predictors,collapse="+")
    
    out5 <- data.frame(marker_name=marker_name,
                       data_name=data_name,
                       type = "Adjusted",
                       model = "Time dependent Cox PH",
                       HR = hr,
                       LB = lb,
                       UB = ub,
                       PVAL = pval,
                       PRED = pred)
    
    
    #########################
    ### Adjusted Simple Cox
    #########################
    hr <- exp(as.numeric(fit$simplecox_model$coefficients["marker1"]))
    wher <- row.names(confint(fit$simplecox_model)) %in% "marker1"
    lb <- exp(as.numeric(confint(fit$simplecox_model)[wher,][1]))
    ub <- exp(as.numeric(confint(fit$simplecox_model)[wher,][2]))
    hr <- round(hr,2)
    lb <- round(lb,2)
    ub <- round(ub,2)
    wher <- row.names(summary(fit$simplecox_model)$coefficients) %in% "marker1"
    pval <- summary(fit$simplecox_model)$coefficients[wher,"Pr(>|z|)"]
    pval <- round(pval,3)
    pval <- format(pval,nsmall=3)
    pred <- paste(fit$simplecox_model_predictors,collapse="+")
    
    out6 <- data.frame(marker_name=marker_name,
                       data_name=data_name,
                       type = "Adjusted",
                       model = "Simple Cox PH",
                       HR = hr,
                       LB = lb,
                       UB = ub,
                       PVAL = pval,
                       PRED = pred)
    
    out <- rbind(out1,out2,out3,out4,out5,out6)
    out
  })
  
  OUT <- do.call("rbind",temp)
  
  if (!is.null(saveopt)) {
    save(OUT,file = paste(dir_out,saveopt,".RData",sep="") )
    write.csv(OUT,file = paste(dir_out,,saveopt,".csv",sep=""),row.names = FALSE)
  }

  OUT
}


