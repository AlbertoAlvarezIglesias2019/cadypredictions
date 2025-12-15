#' @title Predicts event risk using Partly Conditional, Time-Dependent, and Simple Baseline Cox Models
#'
#' @description This function loads longitudinal data, joins it with a separate training/test population set,
#' fits three distinct survival models—a Partly Conditional Cox (PC.Cox) model, a traditional Time-Dependent Cox
#' (TD-Cox) model, and a simple Cox model using only baseline marker value—and generates comparative risk predictions
#' for the test set population.
#' 
#' A key feature is the dynamic reduction of covariates (`Predictors`) in the model fitting process whenever a 
#' convergence warning (often due to perfect separation or collinearity) is detected during the Cox model fitting.
#' 
#' @param datos An optional data frame containing the primary longitudinal data.
#' If provided (\code{NULL} is the default), this data is used instead of loading
#' the file specified by \code{data_name}.
#' @param data_name A character string specifying the prefix of the primary longitudinal data file name.
#' This parameter is only used if \code{datos} is \code{NULL}. The function expects the file at the hardcoded path:
#' "M:/CRF/ICORG/Studies/CADY/Clinical_Study_Report/Report/data/".
#' @param marker_name A character vector string specifying the name of the longitudinal marker(s) column (e.g., "NT_pro_BNP")
#' used as a time-dependent predictor (there could be more than one biomarker). 
#' @param Predictors A character vector specifying names of time-independent covariates (e.g., "Age", "treatment_reg")
#' to be included in all survival models. The list of predictors is dynamically reduced if convergence issues are encountered.
#' @param log_marker A logical value. If \code{TRUE}, the marker variable is transformed using \code{log2()} before
#' model fitting and interpolation. Defaults to \code{FALSE}.
#' @param pred_from An integer specifying the time point (e.g., days) at which the marker value is interpolated,
#' serving as the start time for the risk prediction. Defaults to 150.
#' @param pred_to An integer specifying the total prediction time (e.g., days). The prediction window is
#' \code{pred_to} - \code{pred_from}. Defaults to 240.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{pred_data}: A data frame containing the SubjectID, selected predictors,
#'         event status, the prediction time window (\code{predict_from} to \code{predict_to}),
#'         and the calculated risks from the three models: \code{Risk_PC}, \code{Risk_TDcox}, and \code{Risk_cox_simple}.
#'   \item \code{pc_model}, \code{tdcox_model}, \code{simplecox_model}: The fitted model objects
#'         for the Partly Conditional Cox, Time-Dependent Covariate Cox, and Simple Baseline Cox models, respectively.
#'   \item \code{pc_model_predictors}, \code{tdcox_model_predictors}, \code{simplecox_model_predictors}:
#'         The final character vectors of predictors used in each model after automatic reduction, if any.
#' }
#'
#' @details
#' This function performs a rigorous comparison of three survival prediction methods.
#'
#' \strong{Data Handling and Training/Test Split:}
#' \itemize{
#' \item Data can be provided directly via the \code{datos} argument or loaded from the hardcoded path (via \code{data_name}).
#' Subject assignment to "training" or "test" sets is always read from a separate file ("population_set.csv") at the hardcoded path.
#' \item All three models are trained exclusively on the "training" set.
#' \item Risk predictions are generated exclusively for subjects in the "test" set.
#' }
#'
#' \strong{Marker Interpolation at \code{pred_from} (for PC.Cox and TD-Cox):}
#' Marker values for all subjects are interpolated or carried forward to the prediction start time (\code{pred_from}).
#' \itemize{
#' \item \strong{Interpolation:} Used if measurements exist immediately before and after \code{pred_from} (linear interpolation).
#' \item \strong{Carry Forward:} Used if only measurements before \code{pred_from} are available (the last value is used).
#' \item \strong{Exclusion:} Subjects whose first available marker sample is \emph{after} \code{pred_from} are excluded from the dynamic prediction cohorts.
#' }
#'
#' \strong{Model Details:}
#' \itemize{
#' \item \strong{Partly Conditional Cox (PC.Cox):} Trained on the full longitudinal marker data and uses
#' \code{log2(time_to_sample+1)} as a predictor to model the time dependency of the marker effect.
#' \item \strong{Time-Dependent Cox (TD-Cox):} Uses the \code{survival::tmerge} approach to create time-dependent
#' data, where the marker value is treated as constant between measurements. The marker value interpolated at
#' \code{pred_from} is added to the data for prediction purposes.
#' \item \strong{Simple Baseline Cox:} This model is fitted using only the **baseline** value of the specified marker
#' (assumed to be stored in the column \code{marker_name}_bl). Prediction for this model is the absolute risk of event
#' between \code{pred_from} and \code{pred_to} (i.e., $\hat{S}(\text{pred\_from}) - \hat{S}(\text{pred\_to})$).
#' }
#'
#' @importFrom dplyr filter mutate group_by arrange slice_tail slice ungroup select distinct left_join full_join if_else
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
#' data_name = "cady_data_mp_50",
#' marker_name = "NT_pro_BNP",
#' Predictors = c("Age", "treatment_reg"),
#' log_marker = TRUE,
#' pred_from = 180, # Interpolate marker at 180 days
#' pred_to = 365 # Predict risk up to 365 days
#' )
#' print(results_comparison)
#'
#' # Example access to the prediction data frame:
#' head(results_comparison$pred_data)
#' }
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

  if (!("set" %in% names(dat))) {
    po <- "M:/CRF/ICORG/Studies/CADY/Clinical_Study_Report/Report/data/population_set.csv"
    dat1 <- read.csv(po)
    dat1 <- dat1 %>% select(SubjectID = Label,set)
    dat <- left_join(dat,dat1,by="SubjectID")
  }

  
  DD <- dat
  temp <- as.data.frame(DD[,marker_name])
  names(temp) <- paste("marker",1:length(marker_name),sep="")
  DD <- cbind(DD,temp)
  marker_name_temp <- names(temp)
  #DD$marker <- DD[,c("SubjectID",marker_name)]


  vars <- c("SubjectID",Predictors,marker_name_temp,"time_to_event","status","time_to_sample","set")
  DD <- DD[,vars]


  for (ii in marker_name_temp) DD <- DD[!is.na(DD[,ii]),]
  DD <- DD %>% filter(!is.na(time_to_event))

  if (log_marker) {
    for (ii in marker_name_temp) {
      tt <- min(DD[,ii][DD[,ii]>0],na.rm=TRUE)
      DD[,ii] <- if_else(DD[,ii]==0,log2(DD[,ii]+tt/2),log2(DD[,ii]))
      #DD <- DD %>% mutate(marker = if_else(marker==0,log2(marker+tt/2),log2(marker))  )
    }

  } 

  masterD <- DD
  
  ########################################
  #### Fit a partly conditional Cox model
  ########################################
  DD <- masterD %>% 
    filter(time_to_sample>=0) %>% mutate(log_time_to_sample = log2(time_to_sample+1)) %>% 
    filter(time_to_event>=time_to_sample)%>% 
    filter(set=="training")
  #pccox <-  PC.Cox(
  #  id = "SubjectID",
  #  stime = "time_to_event",
  #  status = "status",
  #  measurement.time = "time_to_sample",  ##survival time and measurement time must be on the same scale.
  #  predictors =c("log_time_to_sample", marker_name_temp,Predictors),
  #  data = DD)
  
  ### CHECK if it gets a warning and reduce the number of predictors
  reploop <- TRUE
  temppred <- Predictors
  while (reploop) {
    pccox <- tryCatch({
      PC.Cox(
        id = "SubjectID",
        stime = "time_to_event",
        status = "status",
        measurement.time = "time_to_sample",  ##survival time and measurement time must be on the same scale.
        predictors =c("log_time_to_sample", marker_name_temp,temppred),
        data = DD)
    }, warning = function(w) {
      TRUE
    })
    if (class(pccox)=="logical") {
      temppred <- temppred[-length(temppred)]
      reploop<-TRUE
    } else {
      reploop<-FALSE}
  }


  fit_pccox <- pccox$model.fit
  fit_pccox_pred <- temppred 
  
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
  D1 <- dd %>% pivot_longer(contains("marker"))
  dd1 <- D1 %>% filter(n>1) %>% group_by(SubjectID,name) %>% mutate(marker_pf = inter(time_to_sample,value))
  # Carry forward if they are only before
  dd2 <- D1 %>% filter(n==1) %>% mutate(marker_pf = value)
  dd <- rbind(dd1,dd2) %>% select(-value)
  dd <- dd %>%
    mutate(time_to_sample=pred_from,
           log_time_to_sample = log2(pred_from),
           marker = marker_pf) 

  nd <- dd %>% select(-when,-n,-marker_pf) %>% distinct() %>% 
    pivot_wider(names_from = name,values_from = marker)
  


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

  out$marker_name <- paste(marker_name,collapse="; ")
  vars <- c("SubjectID",Predictors,marker_name_temp,"marker_name","time_to_event","status","time_to_sample","Risk_PC")

  out_risk_pc <- out[,vars]
  
  
  
  #***************************
  #***
  #*** COX model TD Covariates
  #***
  #***************************
  #if (is.null(datos)) {
  #  pa <- paste("M:/CRF/ICORG/Studies/CADY/Clinical_Study_Report/Report/data/",data_name,".csv",sep="")
  #  dat <- read.csv(pa)
  #} else {dat <- datos}
  
  #if (!("set" %in% names(dat))) {
  #  po <- "M:/CRF/ICORG/Studies/CADY/Clinical_Study_Report/Report/data/population_set.csv"
  #  dat1 <- read.csv(po)
  #  dat1 <- dat1 %>% select(SubjectID = Label,set)
  #  dat <- left_join(dat,dat1,by="SubjectID")
  #}
  
  
  
  #DD <- dat
  #DD$marker <- DD[,marker_name]
  
  
  #vars <- c("SubjectID",Predictors,"marker","time_to_event","status","time_to_sample","set")
  #DD <- DD[,vars]
  
  
  #DD <- DD %>% filter(!is.na(marker))
  #DD <- DD %>% filter(!is.na(time_to_event))
  
  #if (log_marker) {
  #  tt <- min(DD$marker[DD$marker>0],na.rm=TRUE)
  #  DD <- DD %>% mutate(marker = if_else(marker==0,log2(marker+tt/2),log2(marker))  )
  #} 
  
  #masterD <- DD
  
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
  # interpolate if they have before and after
  x <- c(84,168)
  y <- c(130.40,128.70)
  inter <- function(x,y) {
    sl <- (y[2] - y[1])/(x[2] - x[1])
    y[1] + sl*(pred_from-x[1])
  }
  #inter(x,y)
  D1 <- dd %>% pivot_longer(contains("marker"))
  dd1 <- D1 %>% filter(n>1) %>% group_by(SubjectID,name) %>% mutate(marker_pf = inter(time_to_sample,value))
  # Carry forward if they are only before
  dd2 <- D1 %>% filter(n==1) %>% mutate(marker_pf = value)
  dd <- rbind(dd1,dd2) %>% select(-value)
  dd <- dd %>%
    mutate(time_to_sample=pred_from,
           #log_time_to_sample = log2(pred_from),
           marker = marker_pf) 
  
  dd <- dd %>% select(-when,-n,-marker_pf) %>% distinct() %>% 
    pivot_wider(names_from = name,values_from = marker)
  
  masterD <- rbind(masterD,dd) %>% arrange(SubjectID,time_to_sample)
  
  
  ##############################################
  ###Set up data for time dependent Cox model
  ##############################################
  dat1 <- masterD[,c("SubjectID",Predictors,"time_to_event", "status","set")] %>% distinct()
  dat2 <- masterD[,c("SubjectID",marker_name_temp,"time_to_sample")] 
  
  ddd1 <- tmerge(dat1,dat1,id=SubjectID,death=event(time_to_event,status))
  
  texto <- paste(paste(marker_name_temp,"=tdc(time_to_sample,",marker_name_temp,")",sep=""),collapse=",")
  texto <- paste("masterD <- tmerge(ddd1,dat2,id = SubjectID,",texto,")",sep="")
  masterD <- eval(parse(text = texto))
  #masterD <- tmerge(ddd1,dat2,id = SubjectID,marker = tdc(time_to_sample,marker)) 
  
  ddd <- masterD %>% filter(set=="training")
  
  #### OLD
  
  #if (is.null(Predictors)) {
  #  formu <- as.formula(paste("Surv(tstart,tstop,death)~",paste(marker_name_temp,collapse="+"),sep="")) } else {
  #    formu <- as.formula(paste(paste("Surv(tstart,tstop,death)~",paste(marker_name_temp,collapse="+"),sep=""),"+",paste(Predictors,collapse="+"),sep=""))
  #  }
  #fit_tdcox <- coxph(formu,data=ddd)
  
  #### OLD
  
  ### CHECK if it gets a warning and reduce the number of predictors
  reploop <- TRUE
  temppred <- Predictors
  while (reploop) {
    fit_tdcox <- tryCatch({
      
      if (is.null(Predictors)) {
        formu <- as.formula(paste("Surv(tstart,tstop,death)~",paste(marker_name_temp,collapse="+"),sep="")) } else {
          formu <- as.formula(paste(paste("Surv(tstart,tstop,death)~",paste(marker_name_temp,collapse="+"),sep=""),"+",paste(temppred,collapse="+"),sep=""))
        }
      
      coxph(formu,data=ddd)
      
    }, warning = function(w) {
      TRUE
    })
    if (class(fit_tdcox)=="logical") {
      temppred <- temppred[-length(temppred)]
      reploop<-TRUE
    } else {
      reploop<-FALSE}
  }
  
  fit_tdcox_pred <- temppred 
  
  
  
  ###################################################
  ### Make the predictions from pred_from to pred_to
  ###################################################
  ndd <- masterD %>% 
    filter(set=="test")%>%
    filter(tstart==pred_from) 
  
  aa <- summary(survfit(fit_tdcox,newdata = ndd,se.fit=FALSE),times = pred_to )
  ndd$Risk_TDcox <- as.numeric(1-aa$surv)
  out_TDcox <- ndd %>% select(SubjectID,Risk_TDcox)
  
  
  
  
  #******************************
  #***
  #*** COX model SIMPLE baseline
  #***
  #******************************
  
  if (is.null(datos)) {
    pa <- paste("M:/CRF/ICORG/Studies/CADY/Clinical_Study_Report/Report/data/",data_name,".csv",sep="")
    dat <- read.csv(pa)
  } else {dat <- datos}
  
  if (!("set" %in% names(dat))) {
    po <- "M:/CRF/ICORG/Studies/CADY/Clinical_Study_Report/Report/data/population_set.csv"
    dat1 <- read.csv(po)
    dat1 <- dat1 %>% select(SubjectID = Label,set)
    dat <- left_join(dat,dat1,by="SubjectID")
  }
  
  
  
  DD <- dat
  #DD$marker <- DD[,paste(marker_name,"_bl",sep="") ]
  for (ii in 1:length(marker_name_temp) ) {
    DD[,marker_name_temp[ii] ] <- DD[,paste(marker_name[ii],"_bl",sep="") ]
  }
  
  vars <- c("SubjectID",Predictors,marker_name_temp,"time_to_event","status","set")
  DD <- DD[,vars] %>% unique()
  
  
  
  #DD <- DD %>% filter(!is.na(marker))
  for (ii in marker_name_temp) DD <- DD[!is.na(DD[,ii]),]
  DD <- DD %>% filter(!is.na(time_to_event))
  
  #if (log_marker) {
  #  tt <- min(DD$marker[DD$marker>0],na.rm=TRUE)
  #  DD <- DD %>% mutate(marker = if_else(marker==0,log2(marker+tt/2),log2(marker))  )
  #} 
  if (log_marker) {
    for (ii in marker_name_temp) {
      tt <- min(DD[,ii][DD[,ii]>0],na.rm=TRUE)
      DD[,ii] <- if_else(DD[,ii]==0,log2(DD[,ii]+tt/2),log2(DD[,ii]))
      #DD <- DD %>% mutate(marker = if_else(marker==0,log2(marker+tt/2),log2(marker))  )
    }
    
  } 
  
  
  masterD <- DD
  
  ddd <- masterD %>% filter(set=="training")
  
  
  #### OLD
  
  #if (is.null(Predictors)) {
  #  formu <- as.formula(paste("Surv(time_to_event,status)~",paste(marker_name_temp,collapse="+"),sep="")) } else {
  #    formu <- as.formula(paste("Surv(time_to_event,status)~",paste(marker_name_temp,collapse="+"),"+",paste(Predictors,collapse="+"),sep=""))
  #  }
  
  #fit_cox_simple <- coxph(formu,data=ddd)
  
  #### OLD
  
  ### CHECK if it gets a warning and reduce the number of predictors
  reploop <- TRUE
  temppred <- Predictors
  while (reploop) {
    fit_cox_simple <- tryCatch({
      
      if (is.null(Predictors)) {
        formu <- as.formula(paste("Surv(time_to_event,status)~",paste(marker_name_temp,collapse="+"),sep="")) } else {
          formu <- as.formula(paste(paste("Surv(time_to_event,status)~",paste(marker_name_temp,collapse="+"),sep=""),"+",paste(temppred,collapse="+"),sep=""))
        }
      
      coxph(formu,data=ddd)
      
    }, warning = function(w) {
      TRUE
    })
    if (class(fit_cox_simple)=="logical") {
      temppred <- temppred[-length(temppred)]
      reploop<-TRUE
    } else {
      reploop<-FALSE}
  }
  
  fit_cox_simple_pred <- temppred 
  
  
  
  ###################################################
  ### Make the predictions from pred_from to pred_to
  ###################################################
  ndd <- masterD %>% 
    filter(set=="test")
  
  #aa <- summary(survfit(fit_cox_simple,newdata = ndd,se.fit=FALSE),times = pred_to )
  aa1 <- summary(survfit(fit_cox_simple,newdata = ndd,se.fit=FALSE),times = pred_from)
  aa2 <- summary(survfit(fit_cox_simple,newdata = ndd,se.fit=FALSE),times = pred_to)
  ndd$Risk_cox_simple <- as.numeric(aa1$surv-aa2$surv)
  out_simpleCox <- ndd %>% select(SubjectID,Risk_cox_simple)
  
  
  
  
  
  
  

  out <- out_risk_pc %>% full_join(out_TDcox,by="SubjectID") %>% full_join(out_simpleCox,by="SubjectID")
  
  varia <- c("marker_name","SubjectID",Predictors,"time_to_event","status","time_to_sample",marker_name_temp,"Risk_PC","Risk_TDcox","Risk_cox_simple")
  
  out <- out[,varia]
  
  out <- out %>% mutate(predict_to = pred_to)
  out <- out %>% select(marker_name:status,predict_from = time_to_sample,predict_to,all_of(marker_name_temp),Risk_PC,Risk_TDcox,Risk_cox_simple)
  
  
  list(pred_data = out,
       pc_model = fit_pccox,pc_model_predictors = fit_pccox_pred,
       tdcox_model = fit_tdcox,tdcox_model_predictors = fit_tdcox_pred,
       simplecox_model = fit_cox_simple,simplecox_model_predictors = fit_cox_simple_pred)
}

