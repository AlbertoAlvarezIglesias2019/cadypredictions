#' @title Predicts event risk using a Partly Conditional Cox Model with marker interpolation
#'
#' @description This function loads longitudinal patient data, fits a Partly Conditional Cox (PC.Cox) model,
#' interpolates a longitudinal marker value at a specific prediction start time (\code{pred_from}), and then
#' uses the fitted model to predict the event risk at a future time (\code{pred_to}).
#'
#' @param data_name A character string specifying the prefix of the data file name. The function assumes
#'   the file is located at the hardcoded path "M:/CRF/ICORG/Studies/CADY/Clinical_Study_Report/Report/data/".
#' @param marker_name A character string specifying the name of the longitudinal marker column (e.g., "NT_pro_BNP")
#'   in the data frame that should be used as a time-dependent predictor.
#' @param Predictors A character vector specifying names of time-independent covariates (e.g., "Age", "treatment_reg")
#'   to be included in the PC.Cox model.
#' @param log_marker A logical value. If \code{TRUE}, the marker variable is transformed using \code{log2()} before
#'   model fitting and interpolation. Defaults to \code{FALSE}.
#' @param pred_from An integer specifying the time point (e.g., days) at which the marker value should be
#'   interpolated and the risk prediction should begin. Defaults to 150.
#' @param pred_to An integer specifying the time point (e.g., days) to which the risk is predicted. The
#'   prediction window is \code{pred_to} - \code{pred_from}. Defaults to 240.
#'
#' @return A data frame (tibble) containing the original subject data (\code{SubjectID}, \code{time_to_event},
#'   covariates, etc.) along with the predicted risk, named "Risk_PC".
#'
#' @details
#' The function relies on data located at a hardcoded path: "M:/CRF/ICORG/Studies/CADY/Clinical_Study_Report/Report/data/".
#' This path must be accessible for the function to run.
#'
#' \strong{Marker Interpolation:}
#' The marker value at \code{pred_from} is calculated as follows:
#' \itemize{
#'   \item \strong{Interpolation:} For subjects with measurements immediately before and after \code{pred_from}, the marker value
#'     is linearly interpolated.
#'   \item \strong{Carry Forward:} For subjects with measurements only before \code{pred_from}, the last measured value is carried
#'     forward to \code{pred_from}.
#'   \item \strong{Exclusion:} Subjects whose first available marker sample is \emph{after} \code{pred_from} are excluded from the analysis.
#' }
#'
#' The PC.Cox model uses \code{log2(time_to_sample+1)} as a time-dependent predictor in addition to the marker and specified covariates.
#'
#'
#' @examples
#' \dontrun{
#' # Assuming 'cady_data_mp_50.csv' is accessible at the hardcoded path:
#' # Run with default settings (NT_pro_BNP, Age, treatment_reg)
#' result <- predict_risk(
#'   data_name = "cady_data_mp_50",
#'   marker_name = "NT_pro_BNP",
#'   Predictors = c("Age", "treatment_reg"),
#'   log_marker = TRUE,
#'   pred_from = 150,
#'   pred_to = 365
#' )
#' print(result)
#' }
#' @export

predict_risk <- function(data_name = "cady_data_mp_50",
                         marker_name = "NT_pro_BNP",
                         Predictors = c("Age","treatment_reg"),
                         log_marker = FALSE,
                         pred_from = 150,
                         pred_to = 240) {

  pa <- paste("M:/CRF/ICORG/Studies/CADY/Clinical_Study_Report/Report/data/",data_name,".csv",sep="")
  dat <- read.csv(pa)

  DD <- dat
  DD$marker <- DD[,marker_name]


  vars <- c("SubjectID",Predictors,"marker","time_to_event","status","time_to_sample")
  DD <- DD[,vars]


  DD <- DD %>% filter(!is.na(marker))
  DD <- DD %>% filter(!is.na(time_to_event))

  if (log_marker) DD <- DD %>% mutate(marker = log2(marker))

  ########################################
  #### Fit a partly conditional Cox model
  ########################################
  DD <- DD %>% filter(time_to_sample>=0) %>% mutate(log_time_to_sample = log2(time_to_sample+1))
  pccox <-  PC.Cox(
    id = "SubjectID",
    stime = "time_to_event",
    status = "status",
    measurement.time = "time_to_sample",  ##survival time and measurement time must be on the same scale.
    predictors =c("log_time_to_sample", "marker",Predictors),
    data = DD)

  ################
  #**
  ### PREDICTIONS
  #**
  ################

  #########################################################
  ## Pick only biomarker values before and after pred_from
  #########################################################
  dd1 <- DD %>% filter(time_to_sample<= pred_from) %>%
    group_by(SubjectID) %>%
    arrange(time_to_sample) %>%
    slice_tail(n=1) %>%
    ungroup() %>%
    mutate(when = "before")
  dd2 <- DD %>% filter(time_to_sample> pred_from) %>%
    group_by(SubjectID) %>%
    arrange(time_to_sample) %>%
    slice(1) %>%
    ungroup()%>%
    mutate(when = "after")
  dd <- rbind(dd1,dd2) %>% arrange(SubjectID,time_to_sample)

  ###########################################################
  ### Check those with only values before or after pred_from
  ###########################################################
  dd <- dd %>% group_by(SubjectID) %>% mutate(n=n()) %>% ungroup()

  ##################################################
  ### Remove if biomarker sample is after pred_from
  ##################################################
  dd <- dd %>% filter(!(n==1 & when=="after"))

  ##########################################
  ### Interpolate marker level at pred_from
  ##########################################
  # interpolate if they have before and after
  x <- c(84,168)
  y <- c(130.40,128.70)
  inter <- function(x,y) {
    sl <- (y[2] - y[1])/(x[2] - x[1])
    y[1] + sl*(pred_from-x[1])
  }
  #inter(x,y)
  dd1 <- dd %>% filter(n>1) %>% group_by(SubjectID) %>% mutate(marker_pf = inter(time_to_sample,marker))
  # Carry forward if they are only before
  dd2 <- dd %>% filter(n==1) %>% mutate(marker_pf = marker)
  dd <- rbind(dd1,dd2)
  dd <- dd %>%
    mutate(time_to_sample=pred_from,
           log_time_to_sample = log2(pred_from),
           marker = marker_pf)

  nd <- dd %>% select(-when,-n,-marker_pf) %>%
    distinct()

  ###################################################
  ### Make the predictions from pred_from to pred_to
  ###################################################
  nd <- nd %>% filter(time_to_sample <= time_to_event)
  oouu <- predict(pccox, newdata = nd, prediction.time = pred_to-pred_from)
  wher <- str_detect(names(oouu),"risk_")
  names(oouu)[wher] <- "Risk_PC"

  out <- oouu %>% dplyr::select(-t.star)

  out$marker_name <- marker_name
  vars <- c("SubjectID",Predictors,"marker","marker_name","time_to_event","status","time_to_sample","Risk_PC")


  out
}

