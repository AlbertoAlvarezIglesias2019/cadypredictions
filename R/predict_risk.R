#' @title Predicts event risk using Partly Conditional and Time-Dependent Cox Models
#'
#' @description This function loads longitudinal data, joins it with a separate training/test population set,
#' fits two distinct survival models—a Partly Conditional Cox (PC.Cox) model and a traditional Time-Dependent Cox
#' (TD-Cox) model—and generates comparative risk predictions for the test set population.
#'
#' @param data An optional data frame containing the primary longitudinal data.
#'   If provided (\code{NULL} is the default), this data is used instead of loading
#'   the file specified by \code{data_name}.
#' @param data_name A character string specifying the prefix of the primary longitudinal data file name.
#'   This parameter is only used if \code{data} is \code{NULL}. The function expects the file at the hardcoded path:
#'   "M:/CRF/ICORG/Studies/CADY/Clinical_Study_Report/Report/data/".
#' @param marker_name A character string specifying the name of the longitudinal marker column (e.g., "NT_pro_BNP")
#'   used as a time-dependent predictor.
#' @param Predictors A character vector specifying names of time-independent covariates (e.g., "Age", "treatment_reg")
#'   to be included in both survival models.
#' @param log_marker A logical value. If \code{TRUE}, the marker variable is transformed using \code{log2()} before
#'   model fitting and interpolation. Defaults to \code{FALSE}.
#' @param pred_from An integer specifying the time point (e.g., days) at which the marker value is interpolated,
#'   serving as the start time for the risk prediction. Defaults to 150.
#' @param pred_to An integer specifying the total prediction time (e.g., days). The prediction window is
#'   \code{pred_to} - \code{pred_from}. Defaults to 240.
#'
#' @return A list containing:
#' \item{pred_data}{A data frame (tibble) containing subject-level data, the interpolated marker value at \code{pred_from},
#'   and two distinct risk predictions for the test set: \code{Risk_PC} (from PC.Cox) and \code{Risk_TDcox} (from TD-Cox).}
#' \item{pc_model}{The fitted Partly Conditional Cox model object on the training set.}
#' \item{tdcox_model}{The fitted Time-Dependent Cox model object on the training set.}
#' 
#' @details
#' This function performs a rigorous comparison of two dynamic prediction methods.
#'
#' \strong{Data Handling and Training/Test Split:}
#' \itemize{
#'   \item Data can be provided directly via the \code{data} argument or loaded from the hardcoded path (via \code{data_name}).
#'     Subject assignment to "training" or "test" sets is always read from a separate file ("population_set.csv") at the hardcoded path.
#'   \item Both the PC.Cox and TD-Cox models are trained exclusively on the "training" set.
#'   \item Risk predictions are generated exclusively for subjects in the "test" set.
#' }
#'
#' \strong{Marker Interpolation at \code{pred_from}:}
#' Marker values for all subjects are interpolated or carried forward to the prediction start time (\code{pred_from}).
#' \itemize{
#'   \item \strong{Interpolation:} Used if measurements exist immediately before and after \code{pred_from}.
#'   \item \strong{Carry Forward:} Used if only measurements before \code{pred_from} are available (the last value is used).
#'   \item \strong{Exclusion:} Subjects whose first available marker sample is \emph{after} \code{pred_from} are excluded.
#' }
#'
#' \strong{Model Details:}
#' \itemize{
#'   \item \strong{Partly Conditional Cox (PC.Cox):} Trained on the full longitudinal marker data and uses
#'     \code{log2(time_to_sample+1)} as a predictor.
#'   \item \strong{Time-Dependent Cox (TD-Cox):} Uses the \code{survival::tmerge} approach to create time-dependent
#'     data, where the marker value is treated as constant between measurements. The marker value interpolated at
#'     \code{pred_from} is added to the data for prediction purposes.

#' }
#'
#' @importFrom dplyr filter mutate group_by arrange slice_tail slice ungroup select distinct left_join full_join
#' @importFrom stringr str_detect
#' @importFrom partlyconditional PC.Cox
#' @importFrom survival coxph Surv survfit tmerge tdc
#'
#' @examples
#' \dontrun{
#' # Note: This example requires access to the specified hardcoded file paths.
#' # The model is trained on the 'training' set and predictions are made on the 'test' set.
#'
#' results_comparison <- predict_risk(
#'   data_name = "cady_data_mp_50",
#'   marker_name = "NT_pro_BNP",
#'   Predictors = c("Age", "treatment_reg"),
#'   log_marker = TRUE,
#'   pred_from = 180, # Interpolate marker at 180 days
#'   pred_to = 365    # Predict risk up to 365 days
#' )
#' print(results_comparison)
#' 
#' # Note: Requires access to the M: drive path specified internally in the function.
#' # Use a dummy data frame with the required columns for a runnable example.
#'
#' # Create a minimal dummy data set for illustration (cannot run without the M: drive data)
#' dummy_data <- data.frame(
#'   Label = 1:6,
#'   NT_pro_BNP = runif(12, 100, 500),
#'   Age = rep(c(65, 70, 75), each = 2),
#'   treatment_reg = rep(c(1, 0), 3),
#'   time_to_event = rep(c(300, 250, 400), each = 2),
#'   status = rep(c(1, 0, 1), each = 2),
#'   time_to_sample = c(100, 200, 120, 180, 140, 220)
#' )
#' 
#' # You would need a population_set.csv locally for this to run
#' dummy_pop <- data.frame(Label = 1:6, set = c("training", "training", "training", "test", "test", "test"))
#' 
#' # Example call (requires data):
#' result <- predict_risk(
#'   data_name = "your_data_name",
#'   marker_name = "NT_pro_BNP",
#'   Predictors = c("Age", "treatment_reg"),
#'   pred_from = 150,
#'   pred_to = 365
#' )
#' head(result$pred_data)}
#' 
#' @export
#' 
#' 

predict_risk <- function(datos = NULL,
                         data_name = "cady_data_mp_50",
                         marker_name = "NT_pro_BNP",
                         Predictors = NULL,
                         log_marker = FALSE,
                         pred_from = 150,
                         pred_to = 240) {

  if (is.null(datos)) {
    pa <- paste("M:/CRF/ICORG/Studies/CADY/Clinical_Study_Report/Report/data/",data_name,".csv",sep="")
    dat <- read.csv(pa)
  } else {dat <- datos}

  po <- "M:/CRF/ICORG/Studies/CADY/Clinical_Study_Report/Report/data/population_set.csv"
  dat1 <- read.csv(po)
  dat1 <- dat1 %>% select(SubjectID = Label,set)
  dat <- left_join(dat,dat1)
  
  DD <- dat
  DD$marker <- DD[,marker_name]


  vars <- c("SubjectID",Predictors,"marker","time_to_event","status","time_to_sample","set")
  DD <- DD[,vars]


  DD <- DD %>% filter(!is.na(marker))
  DD <- DD %>% filter(!is.na(time_to_event))

  if (log_marker) DD <- DD %>% mutate(marker = log2(marker))

  masterD <- DD
  
  ########################################
  #### Fit a partly conditional Cox model
  ########################################
  DD <- masterD %>% 
    filter(time_to_sample>=0) %>% mutate(log_time_to_sample = log2(time_to_sample+1)) %>% 
    filter(set=="training")
  pccox <-  PC.Cox(
    id = "SubjectID",
    stime = "time_to_event",
    status = "status",
    measurement.time = "time_to_sample",  ##survival time and measurement time must be on the same scale.
    predictors =c("log_time_to_sample", "marker",Predictors),
    data = DD)

  fit_pccox <- pccox$model.fit
  
  ################
  #**
  ### PREDICTIONS
  #**
  ################

  #########################################################
  ## Pick only biomarker values before and after pred_from
  #########################################################
  dd1 <- masterD %>% filter(time_to_sample<= pred_from) %>%
    group_by(SubjectID) %>%
    arrange(time_to_sample) %>%
    slice_tail(n=1) %>%
    ungroup() %>%
    mutate(when = "before")
  dd2 <- masterD %>% filter(time_to_sample> pred_from) %>%
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
  nd <- nd %>%
    filter(time_to_sample <= time_to_event) %>% 
    filter(set=="test")
  
  oouu <- predict(pccox, newdata = nd, prediction.time = pred_to)
  wher <- str_detect(names(oouu),"risk_")
  names(oouu)[wher] <- "Risk_PC"

  out <- oouu %>% dplyr::select(-t.star)

  out$marker_name <- marker_name
  vars <- c("SubjectID",Predictors,"marker","marker_name","time_to_event","status","time_to_sample","Risk_PC")

  out_risk_pc <- out
  
  #***************************
  #***
  #*** COX model
  #***
  #***************************
  
  if (is.null(datos)) {
    pa <- paste("M:/CRF/ICORG/Studies/CADY/Clinical_Study_Report/Report/data/",data_name,".csv",sep="")
    dat <- read.csv(pa)
  } else {dat <- datos}
  
  po <- "M:/CRF/ICORG/Studies/CADY/Clinical_Study_Report/Report/data/population_set.csv"
  dat1 <- read.csv(po)
  dat1 <- dat1 %>% select(SubjectID = Label,set)
  dat <- left_join(dat,dat1)
  
  
  DD <- dat
  DD$marker <- DD[,marker_name]
  
  
  vars <- c("SubjectID",Predictors,"marker","time_to_event","status","time_to_sample","set")
  DD <- DD[,vars]
  
  
  DD <- DD %>% filter(!is.na(marker))
  DD <- DD %>% filter(!is.na(time_to_event))
  
  if (log_marker) DD <- DD %>% mutate(marker = log2(marker))
  
  masterD <- DD
  
  ##############################################
  ### Interpolate the marker value at pred_from
  ##############################################
  dd1 <- masterD %>% filter(time_to_sample<= pred_from) %>%
    group_by(SubjectID) %>%
    arrange(time_to_sample) %>%
    slice_tail(n=1) %>%
    ungroup() %>%
    mutate(when = "before")
  dd2 <- masterD %>% filter(time_to_sample> pred_from) %>%
    group_by(SubjectID) %>%
    arrange(time_to_sample) %>%
    slice(1) %>%
    ungroup()%>%
    mutate(when = "after")
  dd <- rbind(dd1,dd2) %>% arrange(SubjectID,time_to_sample)
  
  ### Check those with only values before or after pred_from
  dd <- dd %>% group_by(SubjectID) %>% mutate(n=n()) %>% ungroup()
  
  ### Remove if biomarker sample is after pred_from
  dd <- dd %>% filter(!(n==1 & when=="after"))
  
  ### interpolate if they have before and after
  dd1 <- dd %>% filter(n>1) %>% group_by(SubjectID) %>% mutate(marker_pf = inter(time_to_sample,marker))
  
  ### Carry forward if they are only before
  dd2 <- dd %>% filter(n==1) %>% mutate(marker_pf = marker)
  
  dd <- rbind(dd1,dd2)
  
  dd <- dd %>%
    mutate(time_to_sample=pred_from,
           marker = marker_pf) %>% 
    select(-when,-n,-marker_pf) %>% unique()
  
  masterD <- rbind(masterD,dd) %>% arrange(SubjectID,time_to_sample)
  
  
  ##############################################
  ###Set up data for time dependent Cox model
  ##############################################
  dat1 <- masterD[,c("SubjectID",Predictors,"time_to_event", "status","set")] %>% distinct()
  dat2 <- masterD[,c("SubjectID","marker","time_to_sample")] 
  
  ddd1 <- tmerge(dat1,dat1,id=SubjectID,death=event(time_to_event,status))
  masterD <- tmerge(ddd1,dat2,id = SubjectID,
                 marker = tdc(time_to_sample,marker)) 
  
  ddd <- masterD %>% filter(set=="training")
  if (is.null(Predictors)) {
    formu <- as.formula("Surv(tstart,tstop,death)~marker")} else {
      formu <- as.formula(paste("Surv(tstart,tstop,death)~marker +",paste(Predictors,collapse="+"),sep=""))
    }
  
  fit_tdcox <- coxph(formu,data=ddd)
  
  
  
  ###################################################
  ### Make the predictions from pred_from to pred_to
  ###################################################
  ndd <- masterD %>% 
    filter(set=="test")%>%
    filter(tstart==pred_from) 
  
  aa <- summary(survfit(fit_tdcox,newdata = ndd,se.fit=FALSE),times = pred_to )
  ndd$Risk_TDcox <- as.numeric(1-aa$surv)
  out <- ndd %>% select(SubjectID,Risk_TDcox)
  

  out <- full_join(out_risk_pc,out)
  
  varia <- c("marker_name","SubjectID",Predictors,"time_to_event","status","time_to_sample","marker","Risk_PC","Risk_TDcox")
  
  out <- out[,varia]
  
  out <- out %>% mutate(predict_to = pred_to)
  out <- out %>% select(marker_name:status,predict_from = time_to_sample,predict_to,marker_at_predfrom = marker,Risk_PC,Risk_TDcox)
  out
  
  list(pred_data = out,pc_model = fit_pccox,tdcox_model = fit_tdcox )
}

