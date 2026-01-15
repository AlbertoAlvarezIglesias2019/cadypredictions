#' Multi-Biomarker Combination Survival Analysis Pipeline
#'
#' @description
#' This function systematically loads various datasets (corresponding to different
#' clinical endpoints), iterates through all possible mathematical combinations of a 
#' user-defined set of biomarkers, and fits six distinct survival models. It estimates 
#' Hazard Ratios (HRs) for each biomarker within the context of its combination.
#'
#' @param dir_in Character string. Path to input CSV data files. 
#'   Expects endpoint files (e.g., \code{dat_nam[i].csv}) and \code{baseline_data.csv}.
#' @param dir_out Character string. Path where \code{RESULTS.RData} and \code{RESULTS.csv} 
#'   will be saved.
#' @param dat_nam Character vector. Base names of data files for different clinical endpoints.
#' @param mar_nam Character vector. The pool of biomarkers from which all 
#'   subsets (combinations) will be generated.
#' @param predictores Character vector. Covariates for "Adjusted" models. 
#'   Must be present in \code{baseline_data.csv}.
#' @param saveopt Character string (optional). File name prefix for saving results.
#' @param ch Logical. If \code{TRUE}, models the "Change from Baseline" (\code{marker - marker_bl}).
#'
#' @return A data frame (\code{OUT}) containing summarized results for every combination 
#'   of data, biomarker subset, adjustment type, and model.
#'
#' @details
#' \strong{Workflow:}
#' \enumerate{
#'   \item **Data Processing:**
#'     \itemize{
#'       \item Loads endpoint data and \code{baseline_data.csv}.
#'       \item Filters for \code{SubjectID}, markers, baseline markers, and survival columns.
#'       \item Merges endpoint data with baseline covariates using a \code{left_join}.
#'     }
#'
#' 
#'
#'   \item **Model Fitting (using \code{predict_risk}):** Six models are fitted per iteration:
#'     \itemize{
#'       \item **Unadjusted Partly Conditional (PC)**: Only log-transformed biomarkers.
#'       \item **Unadjusted Time-Dependent (TD) Cox**: Only log-transformed biomarkers.
#'       \item **Unadjusted Simple Cox PH**: Only baseline biomarker values.
#'       \item **Adjusted Partly Conditional (PC)**: Biomarkers + \code{predictores}.
#'       \item **Adjusted Time-Dependent (TD) Cox**: Biomarkers + \code{predictores}.
#'       \item **Adjusted Simple Cox PH**: Baseline biomarkers + \code{predictores}.
#'     }
#'
#'   \item **Result Extraction:** #'     \itemize{
#'       \item Extracts Hazard Ratios ($\text{HR} = e^{\beta}$), 95\% CIs, and p-values for all 
#'             biomarkers in the current combination (mapped to \code{marker1}, \code{marker2}, etc.).
#'       \item Rounds HR/CIs to 2 decimal places and P-values to 3.
#'     }
#' }
#'
#' 
#'
#' @import cadypredictions
#' @import tidyverse
#' @import partlyconditional
#' @import survival
#'
#' @seealso
#' The core modeling is performed by \code{\link[cadypredictions]{predict_risk}}. See also \code{\link[cadypredictions]{RESULTS}}
#'
#' @examples
#' \dontrun{
#' # Run the analysis for two specific data sets, one marker, and custom predictors
#' custom_predictors <- c("Age", "treatment_reg")
#' RESLIST(
#'   dir_in = "M:/CRF/ICORG/Studies/CADY/Clinical_Study_Report/Report/data/",
#'   dir_out = "M:/CRF/ICORG/Studies/CADY/Clinical_Study_Report/Report/results/",
#'   dat_nam = c("cady_data_mp_50", "cady_data_mp_53"),
#'   mar_nam = c("NT_pro_BNP","BNP),
#'   predictores = custom_predictors
#' )
#' # The results will be saved in "C:/my_results/RESULTS.csv"
#' }
#'
#' @export
#' 
RESLIST <- function(dir_in = "M:/CRF/ICORG/Studies/CADY/Clinical_Study_Report/Report/data/",
                    dir_out = "M:/CRF/ICORG/Studies/CADY/Clinical_Study_Report/Report/results/",
                    dat_nam = c("cady_data_ct","cady_data_drug","cady_data_max_50","cady_data_max_53","cady_data_mp_50","cady_data_mp_53"),
                    mar_nam = c("BNP","NT_pro_BNP","CRP","hsTnI_STAT","Galectin_3"),
                    predictores = c("Age","lvef_mp_bas","diabetes_mellitus_YN","hypertension_YN","dyslipidemia_YN","treatment_reg"),
                    saveopt = NULL,
                    ch = FALSE){
  
  library(cadypredictions)
  library(tidyverse)
  library(partlyconditional)
  library(survival)
  
  # Loop from length 1 to length 4
  all_combs <- lapply(1:length(mar_nam), function(x) combn(mar_nam, x, simplify = FALSE))
  
  # Flatten the list
  all_combs <- unlist(all_combs, recursive = FALSE)
  
  #pointer <- expand.grid(data_name = dat_nam,
  #                       marker_name = mar_nam)
  

  
  datasets_out <- lapply(dat_nam,function (dn){
    
    biomarker_out <- lapply(1:length(all_combs),function(bio){
      
      marker_name <- all_combs[[bio]] 
      marker_name_lab <- marker_name
      if (ch) marker_name_lab <- paste(marker_name," (Change)",sep="")
      
      data_name <- dn
      cat("\n\n\n Row = ",kkk,"; Data and Biomarker: ",data_name," and ",marker_name_lab,"\n\n")  
      
      
      pa <- paste(dir_in,data_name,".csv",sep="")
      masterD <- read.csv(pa)
      #masterD <- masterD %>% select(SubjectID,BNP,NT_pro_BNP,CRP,hsTnI_STAT,Galectin_3,
      #                              BNP_bl,NT_pro_BNP_bl,CRP_bl,hsTnI_STAT_bl,Galectin_3_bl,
      #                              lvef_mp_bas,lvef_max_bas,time_to_event,status,time_to_sample)
      n0 <- "SubjectID"
      n1 <- marker_name
      n2 <- paste(marker_name,"_bl",sep="")
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
      
      nam <- paste("marker",1:length(marker_name),sep="")
      nam_bl <- paste("marker",1:length(marker_name),"_bl",sep="")
      
      wher <- names(fit$pc_model$coefficients) %in% nam
      hr <- exp(as.numeric(fit$pc_model$coefficients[wher]))
      wher <- row.names(confint(fit$pc_model)) %in% nam
      lb <- exp(as.numeric(confint(fit$pc_model)[,1][wher]))
      ub <- exp(as.numeric(confint(fit$pc_model)[,2][wher]))
      hr <- round(hr,2)
      lb <- round(lb,2)
      ub <- round(ub,2)
      wher <- row.names(summary(fit$pc_model)$coefficients) %in% nam
      pval <- summary(fit$pc_model)$coefficients[wher,"Pr(>|z|)"]
      pval <- round(pval,3)
      pval <- format(pval,nsmall=3)
      pred <- paste(fit$pc_model_predictors,collapse="+")
      
      out1 <- data.frame(marker_name=marker_name_lab,
                         data_name=data_name,
                         type = "Unadjusted",
                         model = "Partly conditional",
                         model_num = 1,
                         HR = hr,
                         LB = lb,
                         UB = ub,
                         PVAL = pval,
                         PRED = pred,
                         N = fit$event_counts %>% filter(Model == "pccox",Set=="Training",Counts=="Total sample") %>% pull(Value),
                         E = fit$event_counts %>% filter(Model == "pccox",Set=="Training",Counts=="Total events") %>% pull(Value),
                         EinW = fit$event_counts %>% filter(Model == "pccox",Set=="Training",Counts=="Total events in window") %>% pull(Value))
      
      ###################
      ### Unadjusted TD
      ###################
      wher <- names(fit$tdcox_model$coefficients) %in% nam
      hr <- exp(as.numeric(fit$tdcox_model$coefficients[wher]))
      wher <- row.names(confint(fit$tdcox_model)) %in% nam
      lb <- exp(as.numeric(confint(fit$tdcox_model)[,1][wher]))
      ub <- exp(as.numeric(confint(fit$tdcox_model)[,2][wher]))
      hr <- round(hr,2)
      lb <- round(lb,2)
      ub <- round(ub,2)
      wher <- row.names(summary(fit$tdcox_model)$coefficients) %in% nam
      pval <- summary(fit$tdcox_model)$coefficients[wher,"Pr(>|z|)"]
      pval <- round(pval,3)
      pval <- format(pval,nsmall=3)
      pred <- paste(fit$tdcox_model_predictors,collapse="+")
      
      out2 <- data.frame(marker_name=marker_name_lab,
                         data_name=data_name,
                         type = "Unadjusted",
                         model = "Time dependent Cox PH",
                         model_num = 2,
                         HR = hr,
                         LB = lb,
                         UB = ub,
                         PVAL = pval,
                         PRED = pred,
                         N = fit$event_counts %>% filter(Model == "tdcox",Set=="Training",Counts=="Total sample") %>% pull(Value),
                         E = fit$event_counts %>% filter(Model == "tdcox",Set=="Training",Counts=="Total events") %>% pull(Value),
                         EinW = fit$event_counts %>% filter(Model == "tdcox",Set=="Training",Counts=="Total events in window") %>% pull(Value))
      
      
      ##########################
      ### Unadjusted Simple Cox
      ##########################
      wher <- names(fit$simplecox_model$coefficients) %in% nam_bl
      hr <- exp(as.numeric(fit$simplecox_model$coefficients[wher]))
      wher <- row.names(confint(fit$simplecox_model)) %in% nam_bl
      lb <- exp(as.numeric(confint(fit$simplecox_model)[,1][wher]))
      ub <- exp(as.numeric(confint(fit$simplecox_model)[,2][wher]))
      hr <- round(hr,2)
      lb <- round(lb,2)
      ub <- round(ub,2)
      wher <- row.names(summary(fit$simplecox_model)$coefficients) %in% nam_bl
      pval <- summary(fit$simplecox_model)$coefficients[wher,"Pr(>|z|)"]
      pval <- round(pval,3)
      pval <- format(pval,nsmall=3)
      pred <- paste(fit$simplecox_model_predictors,collapse="+")
      
      out3 <- data.frame(marker_name=marker_name_lab,
                         data_name=data_name,
                         type = "Unadjusted",
                         model = "Simple Cox PH",
                         model_num = 3,
                         HR = hr,
                         LB = lb,
                         UB = ub,
                         PVAL = pval,
                         PRED = pred,
                         N = fit$event_counts %>% filter(Model == "cox_simple",Set=="Training",Counts=="Total sample") %>% pull(Value),
                         E = fit$event_counts %>% filter(Model == "cox_simple",Set=="Training",Counts=="Total events") %>% pull(Value),
                         EinW = fit$event_counts %>% filter(Model == "cox_simple",Set=="Training",Counts=="Total events in window") %>% pull(Value))
      
      
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
      
      wher <- names(fit$pc_model$coefficients) %in% nam
      hr <- exp(as.numeric(fit$pc_model$coefficients[wher]))
      wher <- row.names(confint(fit$pc_model)) %in% nam
      lb <- exp(as.numeric(confint(fit$pc_model)[,1][wher]))
      ub <- exp(as.numeric(confint(fit$pc_model)[,2][wher]))
      hr <- round(hr,2)
      lb <- round(lb,2)
      ub <- round(ub,2)
      wher <- row.names(summary(fit$pc_model)$coefficients) %in% nam
      pval <- summary(fit$pc_model)$coefficients[wher,"Pr(>|z|)"]
      pval <- round(pval,3)
      pval <- format(pval,nsmall=3)
      pred <- paste(fit$pc_model_predictors,collapse="+")
    
      out4 <- data.frame(marker_name=marker_name_lab,
                         data_name=data_name,
                         type = "Adjusted",
                         model = "Partly conditional",
                         model_num = 4,
                         HR = hr,
                         LB = lb,
                         UB = ub,
                         PVAL = pval,
                         PRED = pred,
                         N = fit$event_counts %>% filter(Model == "pccox",Set=="Training",Counts=="Total sample") %>% pull(Value),
                         E = fit$event_counts %>% filter(Model == "pccox",Set=="Training",Counts=="Total events") %>% pull(Value),
                         EinW = fit$event_counts %>% filter(Model == "pccox",Set=="Training",Counts=="Total events in window") %>% pull(Value))
      
      
      
      
      
      ################
      ### Adjusted TD
      ################
      wher <- names(fit$tdcox_model$coefficients) %in% nam
      hr <- exp(as.numeric(fit$tdcox_model$coefficients[wher]))
      wher <- row.names(confint(fit$tdcox_model)) %in% nam
      lb <- exp(as.numeric(confint(fit$tdcox_model)[,1][wher]))
      ub <- exp(as.numeric(confint(fit$tdcox_model)[,2][wher]))
      hr <- round(hr,2)
      lb <- round(lb,2)
      ub <- round(ub,2)
      wher <- row.names(summary(fit$tdcox_model)$coefficients) %in% nam
      pval <- summary(fit$tdcox_model)$coefficients[wher,"Pr(>|z|)"]
      pval <- round(pval,3)
      pval <- format(pval,nsmall=3)
      pred <- paste(fit$tdcox_model_predictors,collapse="+")
      
      out5 <- data.frame(marker_name=marker_name_lab,
                         data_name=data_name,
                         type = "Adjusted",
                         model = "Time dependent Cox PH",
                         model_num = 5,
                         HR = hr,
                         LB = lb,
                         UB = ub,
                         PVAL = pval,
                         PRED = pred,
                         N = fit$event_counts %>% filter(Model == "tdcox",Set=="Training",Counts=="Total sample") %>% pull(Value),
                         E = fit$event_counts %>% filter(Model == "tdcox",Set=="Training",Counts=="Total events") %>% pull(Value),
                         EinW = fit$event_counts %>% filter(Model == "tdcox",Set=="Training",Counts=="Total events in window") %>% pull(Value))
      
      
      #########################
      ### Adjusted Simple Cox
      #########################
      wher <- names(fit$simplecox_model$coefficients) %in% nam_bl
      hr <- exp(as.numeric(fit$simplecox_model$coefficients[wher]))
      wher <- row.names(confint(fit$simplecox_model)) %in% nam_bl
      lb <- exp(as.numeric(confint(fit$simplecox_model)[,1][wher]))
      ub <- exp(as.numeric(confint(fit$simplecox_model)[,2][wher]))
      hr <- round(hr,2)
      lb <- round(lb,2)
      ub <- round(ub,2)
      wher <- row.names(summary(fit$simplecox_model)$coefficients) %in% nam_bl
      pval <- summary(fit$simplecox_model)$coefficients[wher,"Pr(>|z|)"]
      pval <- round(pval,3)
      pval <- format(pval,nsmall=3)
      pred <- paste(fit$simplecox_model_predictors,collapse="+")
      
      out6 <- data.frame(marker_name=marker_name_lab,
                         data_name=data_name,
                         type = "Adjusted",
                         model = "Simple Cox PH",
                         model_num = 6,
                         HR = hr,
                         LB = lb,
                         UB = ub,
                         PVAL = pval,
                         PRED = pred,
                         N = fit$event_counts %>% filter(Model == "cox_simple",Set=="Training",Counts=="Total sample") %>% pull(Value),
                         E = fit$event_counts %>% filter(Model == "cox_simple",Set=="Training",Counts=="Total events") %>% pull(Value),
                         EinW = fit$event_counts %>% filter(Model == "cox_simple",Set=="Training",Counts=="Total events in window") %>% pull(Value))
      
      out <- rbind(out1,out2,out3,out4,out5,out6)
      row.names(out) <- NULL
      mb <- paste(marker_name,collapse=" + ")
      out <- out %>% mutate(model_bio = mb)
      out <- out %>% select(marker_name:model_num,model_bio,HR:EinW)
      
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


