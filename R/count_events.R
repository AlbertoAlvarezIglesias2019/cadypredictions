#' Summarize Sample Sizes and Event Counts Across Data Splits
#'
#' This function calculates descriptive statistics for a survival analysis dataset,
#' partitioned by "Training" and "Test" sets. It provides counts for the total 
#' number of unique subjects, total events, and events occurring within a specific 
#' time window defined by the user.
#'
#' @param MASD A data frame containing the master dataset. It must include the 
#'   following columns: 
#'   \itemize{
#'     \item \code{set}: A character or factor indicating the split ("training" or "test").
#'     \item \code{SubjectID}: Unique identifier for participants.
#'     \item \code{status}: Binary event indicator (1 for event, 0 for censoring).
#'     \item \code{time_to_event}: Numeric time to the event or censoring.
#'   }
#' @param pred_from Numeric. The start of the time window of interest for 
#'   calculating "Events in window".
#' @param pred_to Numeric. The end of the time window of interest for 
#'   calculating "Events in window".
#'
#' @return A data frame with three columns:
#' \itemize{
#'   \item \code{Set}: The data partition ("Overall", "Training", or "Test").
#'   \item \code{Counts}: The metric being measured ("Total sample", "Total events", 
#'     or "Total events in window").
#'   \item \code{Value}: The numeric count for that metric.
#' }
#'
#' @details
#' The function processes the data in three stages:
#' \enumerate{
#'   \item **Training Set:** Filters rows where \code{set == "training"} and calculates counts.
#'   \item **Test Set:** Filters rows where \code{set == "test"} and calculates counts.
#'   \item **Overall:** Aggregates the results from both sets.
#' }
#' Unique subjects are identified using \code{SubjectID} to ensure that total sample 
#' sizes and total event counts are not inflated by longitudinal repetitions of the same subject.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' event_summary <- count_events(MASD = my_data, pred_from = 150, pred_to = 240)
#' print(event_summary)
#' }
#'
#' @import dplyr
#' @export
count_events <- function(MASD,pred_from,pred_to) {
  
  
  #############
  ###
  ### Training 
  ###
  #############
  masd <- MASD %>% filter(set=="training")

  ### Total sample size
  ttt1 <- masd  %>% select(SubjectID) %>% distinct() %>% pull(SubjectID)
  ttt1 <- length(ttt1)
  
  
  ### Total events 
  ttt2 <- masd  %>% select(SubjectID,status) %>% distinct() %>% pull(status)
  ttt2 <- sum(ttt2)
  
  ### Total events in window
  ttt3 <- masd  %>% select(SubjectID,time_to_event,status) %>% distinct() 
  ttt3 <- ttt3 %>% mutate(inwindow = time_to_event>=pred_from & time_to_event<= pred_to)
  ttt3 <- ttt3 %>% filter(inwindow)
  ttt3 <- sum(ttt3$status)
  
  oouutt1 <- data.frame(Set  = "Training",
                        Counts = c("Total sample","Total events","Total events in window"),
                        Value = c(ttt1,ttt2,ttt3))
  
  #############
  ###
  ### Test 
  ###
  #############
  masd <- MASD %>% filter(set=="test")
  
  ### Total sample size
  ttt1 <- masd  %>% select(SubjectID) %>% distinct() %>% pull(SubjectID)
  ttt1 <- length(ttt1)
  
  
  ### Total events 
  ttt2 <- masd  %>% select(SubjectID,status) %>% distinct() %>% pull(status)
  ttt2 <- sum(ttt2)
  
  ### Total events in window
  ttt3 <- masd  %>% select(SubjectID,time_to_event,status) %>% distinct() 
  ttt3 <- ttt3 %>% mutate(inwindow = time_to_event>=pred_from & time_to_event<= pred_to)
  ttt3 <- ttt3 %>% filter(inwindow)
  ttt3 <- sum(ttt3$status)
  
  oouutt2 <- data.frame(Set  = "Test",
                        Counts = c("Total sample","Total events","Total events in window"),
                        Value = c(ttt1,ttt2,ttt3))
  
  
  #############
  ###
  ### Overall 
  ###
  #############
  oouutt <- rbind(oouutt1,oouutt2)
  
  oouutt3 <- oouutt %>% group_by(Counts) %>% summarise(Value = sum(Value))
  oouutt3 <- oouutt3 %>% mutate(Set = "Overall") %>% select(Set,Counts,Value)
  
  oouutt <- rbind(oouutt3,oouutt)
  
  oouutt <- oouutt %>%
    mutate(Counts = ordered(Counts,levels = c("Total sample","Total events","Total events in window"))) %>% 
    mutate(Set = ordered(Set ,levels = c("Overall","Training","Test"))) %>% 
    arrange(Set,Counts)
  
  oouutt
}