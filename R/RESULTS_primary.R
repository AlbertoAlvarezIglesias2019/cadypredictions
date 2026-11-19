#' Run Adjusted Survival Models Across Multiple Biomarkers and Clinical Endpoints
#'
#' This function systematically loads datasets corresponding to different clinical
#' endpoints, iterates through a user-defined set of biomarkers, and fits
#' adjusted survival models using both the Partly Conditional (PC) model and the
#' Time-Dependent (TD) Cox Proportional Hazards model to estimate Hazard Ratios (HRs).
#'
#' @param dir_in A character string specifying the path to the directory where the
#'   input CSV data files are located. Defaults to
#'   \code{"M:/CRF/ICORG/Studies/CADY/Clinical_Study_Report/Report/data/"}. The
#'   function expects files named \code{<dat_nam[i]>.csv} and \code{baseline_data.csv}.
#' @param dir_out A character string specifying the path to the directory where the
#'   final results (\code{<saveopt>.RData} and \code{<saveopt>.csv}) are saved when
#'   \code{saveopt} is provided. Defaults to
#'   \code{"M:/CRF/ICORG/Studies/CADY/Clinical_Study_Report/Report/results/"}.
#' @param dat_nam A character vector listing the base names of the data files
#'   representing different clinical endpoint definitions. Defaults to
#'   \code{c("cady_data_ct")}.
#' @param mar_nam A character vector listing the biomarker names to analyze.
#'   Must correspond to columns in the loaded data (along with \code{<marker>_bl}
#'   baseline columns). Defaults to
#'   \code{c("BNP", "NT_pro_BNP", "CRP", "hsTnI_STAT", "Galectin_3")}.
#' @param predictores A character vector specifying the names of the covariates
#'   to be included in the adjusted models. Any covariates missing from the endpoint
#'   dataset are merged from \code{baseline_data.csv}. Defaults to
#'   \code{c("Age", "lvef_mp_bas", "diabetes_mellitus_YN", "hypertension_YN",
#'   "dyslipidemia_YN", "treatment_reg", "lapatinib_YN", "rt_YN")}.
#' @param saveopt A character string (optional). If provided, the final aggregated
#'   results data frame (\code{OUT}) is saved to \code{dir_out} as both
#'   \code{<saveopt>.RData} and \code{<saveopt>.csv}. If \code{NULL} (default), no
#'   files are written to disk.
#' @param ch Logical. If \code{TRUE}, models the change from baseline
#'   (\code{marker - marker_bl}) and appends \code{" (Change)"} to
#'   \code{marker_name} in the output. Defaults to \code{FALSE}.
#' @param conf Numeric. The confidence level used to compute the custom confidence
#'   interval bounds (\code{LB} and \code{UB}) in addition to the standard 95\%
#'   confidence bounds (\code{LB95} and \code{UB95}). Defaults to \code{0.995}.
#'
#' @return A data frame (\code{OUT}) containing the summarized results for every
#'   combination of dataset, biomarker, and adjusted model, with the following columns:
#' \itemize{
#'   \item \code{marker_name}: The biomarker analyzed (suffixed with \code{" (Change)"} when \code{ch = TRUE}).
#'   \item \code{data_name}: The endpoint dataset name used.
#'   \item \code{type}: Adjustment status (\code{"Adjusted"}).
#'   \item \code{model}: The survival model type (\code{"Partly conditional"} or \code{"Time dependent Cox PH"}).
#'   \item \code{HR}: The estimated Hazard Ratio (\eqn{e^{\beta}}) for the \eqn{\log_2}-transformed biomarker, rounded to 2 decimal places.
#'   \item \code{LB95}, \code{UB95}: The lower and upper bounds of the 95\% Confidence Interval for the HR, rounded to 2 decimal places.
#'   \item \code{LB}, \code{UB}: The lower and upper bounds of the Confidence Interval at level \code{conf} for the HR, rounded to 2 decimal places.
#'   \item \code{PVAL}: The two-sided p-value for the biomarker coefficient, formatted to 3 decimal places.
#'   \item \code{PRED}: A \code{"+"}-delimited string listing the final predictors retained in the fitted model.
#'   \item \code{N}, \code{E}, \code{EinW}: Training-set counts for Total Sample, Total Events, and Total Events in the Prediction Window (\code{150} to \code{240}), respectively.
#' }
#'
#' @details
#' The function creates an expansion grid (\code{pointer}) across all combinations of
#' \code{dat_nam} and \code{mar_nam} and iterates over each pair:
#'
#' \enumerate{
#'   \item \strong{Data Loading & Merging:}
#'     \itemize{
#'       \item Reads \code{<data_name>.csv} from \code{dir_in} and subsets to \code{SubjectID}, biomarker columns, baseline biomarker columns (\code{<marker>_bl}), survival/timing variables (\code{lvef_mp_bas}, \code{lvef_max_bas}, \code{time_to_event}, \code{status}, \code{time_to_sample}), and any available \code{predictores}.
#'       \item Reads \code{baseline_data.csv} and left-joins any requested \code{predictores} not already present in the endpoint file by \code{SubjectID}.
#'     }
#'   \item \strong{Model Fitting (via \code{\link[cadypredictions]{predict_risk}}):}
#'     Calls \code{predict_risk()} with \code{log_marker = TRUE}, \code{pred_from = 150}, \code{pred_to = 240}, and \code{change = ch}, extracting results for two adjusted models:
#'     \itemize{
#'       \item \strong{Adjusted Partly Conditional (PC):} Models the \eqn{\log_2}-transformed biomarker alongside \code{predictores}.
#'       \item \strong{Adjusted Time-Dependent (TD) Cox PH:} Models the \eqn{\log_2}-transformed biomarker alongside \code{predictores}.
#'     }
#'     \emph{(Note: Although \code{out6} is computed internally for the Simple Cox PH model, only the Partly Conditional (\code{out4}) and Time-Dependent Cox PH (\code{out5}) rows are combined into the returned data frame.)}
#'   \item \strong{Result Extraction & Export:}
#'     Extracts the biomarker HR (\code{marker1}), 95\% CIs, custom \code{conf}-level CIs, p-value, final predictor list, and training event counts, binds all rows into \code{OUT}, and optionally writes \code{.RData} and \code{.csv} files to \code{dir_out}.
#' }
#'
#' @import cadypredictions
#' @import tidyverse
#' @import partlyconditional
#' @import survival
#' @importFrom base paste read.csv data.frame round format
#' @importFrom dplyr select left_join filter pull
#' @importFrom stats confint
#'
#' @seealso
#' \code{\link[cadypredictions]{predict_risk}} for the underlying model fitting and event-count extraction.
#'
#' @examples
#' \dontrun{
#' # Run the primary adjusted analysis for two datasets, one marker, and custom predictors
#' custom_predictors <- c("Age", "treatment_reg")
#' res <- RESULTS_primary(
#'   dir_in = "C:/my_data/",
#'   dir_out = "C:/my_results/",
#'   dat_nam = c("cady_data_mp_50", "cady_data_mp_53"),
#'   mar_nam = "NT_pro_BNP",
#'   predictores = custom_predictors,
#'   saveopt = "RESULTS_primary",
#'   conf = 0.995
#' )
#' # Results are returned in `res` and saved to "C:/my_results/RESULTS_primary.csv"
#' }
#'
#' @export
#' 

RESULTS_primary <- function(dir_in = "M:/CRF/ICORG/Studies/CADY/Clinical_Study_Report/Report/data/",
                   dir_out = "M:/CRF/ICORG/Studies/CADY/Clinical_Study_Report/Report/results/",
                   dat_nam = c("cady_data_ct"),
                   mar_nam = c("BNP","NT_pro_BNP","CRP","hsTnI_STAT","Galectin_3"),
                   predictores = c("Age","lvef_mp_bas","diabetes_mellitus_YN","hypertension_YN","dyslipidemia_YN",
                                   "treatment_reg","lapatinib_YN","rt_YN"),
                   saveopt = NULL,
                   ch = FALSE,
                   conf = 0.995){
  
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
    marker_name_lab <- if_else(ch,paste(marker_name," (Change)",sep=""),marker_name)
    
    #data_name <- "cady_data_ct" 
    #marker_name <- "hsTnI_STAT"
    
    cat("\n\n\n Row = ",kkk,"; Data and Biomarker: ",data_name," and ",marker_name_lab,"\n\n")  
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
    
    hr <- exp(as.numeric(fit$pc_model$coefficients["marker1"]))
    wher <- row.names(confint(fit$pc_model)) %in% "marker1"
    lb <- exp(as.numeric(confint(fit$pc_model,level = conf)[wher,][1]))
    ub <- exp(as.numeric(confint(fit$pc_model,level = conf)[wher,][2]))
    lb95 <- exp(as.numeric(confint(fit$pc_model)[wher,][1]))
    ub95 <- exp(as.numeric(confint(fit$pc_model)[wher,][2]))
    hr <- round(hr,2)
    lb95 <- round(lb95,2)
    ub95 <- round(ub95,2)
    lb <- round(lb,2)
    ub <- round(ub,2)
    wher <- row.names(summary(fit$pc_model)$coefficients) %in% "marker1"
    pval <- summary(fit$pc_model)$coefficients[wher,"Pr(>|z|)"]
    pval <- round(pval,3)
    pval <- format(pval,nsmall=3)
    pred <- paste(fit$pc_model_predictors,collapse="+")
    
    out4 <- data.frame(marker_name=marker_name_lab,
                       data_name=data_name,
                       type = "Adjusted",
                       model = "Partly conditional",
                       HR = hr,
                       LB95 = lb95,
                       UB95 = ub95,
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
    hr <- exp(as.numeric(fit$tdcox_model$coefficients["marker1"]))
    wher <- row.names(confint(fit$tdcox_model)) %in% "marker1"
    lb <- exp(as.numeric(confint(fit$tdcox_model,level = conf)[wher,][1]))
    ub <- exp(as.numeric(confint(fit$tdcox_model,level = conf)[wher,][2]))
    lb95 <- exp(as.numeric(confint(fit$pc_model)[wher,][1]))
    ub95 <- exp(as.numeric(confint(fit$pc_model)[wher,][2]))
    hr <- round(hr,2)
    lb95 <- round(lb95,2)
    ub95 <- round(ub95,2)
    lb <- round(lb,2)
    ub <- round(ub,2)
    wher <- row.names(summary(fit$tdcox_model)$coefficients) %in% "marker1"
    pval <- summary(fit$tdcox_model)$coefficients[wher,"Pr(>|z|)"]
    pval <- round(pval,3)
    pval <- format(pval,nsmall=3)
    pred <- paste(fit$tdcox_model_predictors,collapse="+")
    
    out5 <- data.frame(marker_name=marker_name_lab,
                       data_name=data_name,
                       type = "Adjusted",
                       model = "Time dependent Cox PH",
                       HR = hr,
                       LB95 = lb95,
                       UB95 = ub95,
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
    hr <- exp(as.numeric(fit$simplecox_model$coefficients["marker1_bl"]))
    wher <- row.names(confint(fit$simplecox_model)) %in% "marker1_bl"
    lb <- exp(as.numeric(confint(fit$simplecox_model,level = conf)[wher,][1]))
    ub <- exp(as.numeric(confint(fit$simplecox_model,level = conf)[wher,][2]))
    hr <- round(hr,2)
    lb <- round(lb,2)
    ub <- round(ub,2)
    wher <- row.names(summary(fit$simplecox_model)$coefficients) %in% "marker1_bl"
    pval <- summary(fit$simplecox_model)$coefficients[wher,"Pr(>|z|)"]
    pval <- round(pval,3)
    pval <- format(pval,nsmall=3)
    pred <- paste(fit$simplecox_model_predictors,collapse="+")
    
    out6 <- data.frame(marker_name=marker_name_lab,
                       data_name=data_name,
                       type = "Adjusted",
                       model = "Simple Cox PH",
                       HR = hr,
                       LB = lb,
                       UB = ub,
                       PVAL = pval,
                       PRED = pred,
                       N = fit$event_counts %>% filter(Model == "cox_simple",Set=="Training",Counts=="Total sample") %>% pull(Value),
                       E = fit$event_counts %>% filter(Model == "cox_simple",Set=="Training",Counts=="Total events") %>% pull(Value),
                       EinW = fit$event_counts %>% filter(Model == "cox_simple",Set=="Training",Counts=="Total events in window") %>% pull(Value))
    
    out <- rbind(out4,out5)
    out
  })
  
  OUT <- do.call("rbind",temp)
  
  if (!is.null(saveopt)) {
    save(OUT,file = paste(dir_out,saveopt,".RData",sep="") )
    write.csv(OUT,file = paste(dir_out,saveopt,".csv",sep=""),row.names = FALSE)
  }

  OUT
}


