predict.PC_cox <- function(object, newdata,  prediction.time, ...){
  
  stopifnot(is.numeric(prediction.time))
  prediction.time = sort(prediction.time)
  
  stopifnot(is.data.frame(newdata))
  #check if newdata has the right variables
  if(!all(is.element(object$variable.names, names(newdata)) )){
    stop("variable(s): {", paste(object$names[which(!is.element(object$variable.names, names(newdata)))], collapse = ", "), "}, not found in 'newdata' ")
  }
  
  
  #only keep complete cases
  mycomplete <- complete.cases(newdata[,object$variable.names]);
  
  #check for missing data and throw it out, print a warning
  if(nrow(newdata)!=sum(mycomplete)){
    warning(paste(nrow(newdata)-sum(mycomplete), "observation(s) were removed due to missing data \n New sample size is now:", sum(mycomplete)))
    newdata <- newdata[mycomplete,]
    
  }
  
  #remove any observations with survival time less than meas.time
  wonky.times <-  newdata[[object$variable.names[2]]] < newdata[[object$variable.names[4]]]
  if(any(wonky.times)) {
    warning(paste0("... removing ",  sum(wonky.times) , " observations where ", object$variable.names[4], " is greater than outcome variable ", object$variable.names[2], ".\n"))
  }
  
  newdata <- newdata[!wonky.times,]
  
  
  
  my.data <- data.frame(id = newdata[[object$variable.names[1]]],
                        stime = newdata[[object$variable.names[2]]],
                        status = newdata[[object$variable.names[3]]],
                        measurement.time = newdata[[object$variable.names[4]]],
                        t.star = newdata[[object$variable.names[2]]] - newdata[[object$variable.names[4]]])
  
  
  
  names(my.data) <- c(object$variable.names[1:4], "t.star")
  my.data <- cbind(my.data, newdata[,object$variable.names[-c(1:4)], drop = FALSE])
  
  #sort by measurement time
  id <- object$variable.names[[1]]
  measurement.time <-  object$variable.names[[4]]
  my.data <- my.data %>% arrange(id, measurement.time)
  
  
  
  #calculate risk
  # we calculate risk for the test set, last observation before conditioning time
  
  newdata[["t.star"]] <-  newdata[[object$variable.names[2]]] - newdata[[object$variable.names[4]]]
  #return risk
  aa <- summary(survfit( object$model.fit, newdata = newdata, se.fit = F),
                times = prediction.time)
  aa <- as.data.frame(t(1-aa$surv))
  names(aa) <- paste0("risk_", prediction.time)
  
  out <- bind_cols(newdata, aa)
  
  return(out)
  
  
  
  
  
}