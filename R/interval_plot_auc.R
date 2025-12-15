
#' Create an Interval Plot (Forest Plot) for AUC Values
#'
#' This function generates an interval plot, commonly known as a forest plot,
#' to visualize the Area Under the Curve (AUC) and its confidence intervals
#' for a specific biomarker across different study settings or definitions.
#' The plot is created using \code{ggplot2}.
#'
#' @param dat A data frame containing the AUC results. This data frame is
#'   expected to have columns like \code{marker_name}, \code{data_name},
#'   \code{type}, \code{model}, \code{AUC}, \code{LB}, \code{UB},
#'   \code{LB_nobon}, and \code{UB_nobon}.
#' @param marker_nam A character string specifying the single biomarker to plot.
#'   Must match a value in the \code{marker_name} column of \code{dat}.
#' @param boncorrect A logical value. If \code{TRUE} (default), the plot uses
#'   the Bonferroni-corrected confidence intervals (\code{LB} and \code{UB}).
#'   If \code{FALSE}, it uses the uncorrected confidence intervals
#'   (\code{LB_nobon} and \code{UB_nobon}).
#'
#' @return A \code{ggplot} object representing the interval plot.
#'
#' @details
#' The function first processes the input data:
#' \itemize{
#'   \item It creates a human-readable label (\code{lab1}) for the biomarker
#'     from \code{marker_name}.
#'   \item It creates a human-readable label (\code{lab2}) for the CT definition
#'     from \code{data_name}.
#'   \item It combines these into a final y-axis label (\code{lab}) which includes
#'     the CT definition, type, and model.
#'   \item It filters the data to include only the specified \code{marker_nam}.
#'   \item The resulting data is arranged by \code{AUC} and the \code{lab}
#'     column is converted to a factor to ensure the plot is ordered by AUC.
#' }
#' The plot itself is an interval plot showing the AUC estimate as a point and
#' the confidence interval as a horizontal error bar. A vertical dashed line
#' is added at an AUC of 0.5 for reference.
#' 
#'
#' @importFrom dplyr mutate case_when filter arrange
#' @importFrom ggplot2 ggplot aes geom_errorbarh labs geom_point geom_vline theme_minimal theme element_text
#'
#' @examples
#' \dontrun{
#' # Assuming 'auc_data' is a data frame structured as expected
#' # and includes results for 'BNP'
#'
#' # Plot with Bonferroni correction (default)
#' p_bnp_corrected <- interval_plot_auc(
#'   dat = auc_data,
#'   marker_nam = "BNP",
#'   boncorrect = TRUE
#' )
#' print(p_bnp_corrected)
#'
#' # Plot without Bonferroni correction
#' p_bnp_uncorrected <- interval_plot_auc(
#'   dat = auc_data,
#'   marker_nam = "BNP",
#'   boncorrect = FALSE
#' )
#' print(p_bnp_uncorrected)
#' }
#'
#' @export
#' 

interval_plot_auc <- function(dat,marker_nam,boncorrect = TRUE) {
  
  dat <- dat %>% mutate(lab1 = case_when(marker_name=="BNP"~"BNP",
                                         marker_name=="NT_pro_BNP"~"NT-proBNP",
                                         marker_name=="CRP"~"CRP",
                                         marker_name=="hsTnI_STAT"~"hs-Troponin I",
                                         marker_name=="Galectin_3"~"Galectin-3"))
  
  dat <- dat %>% mutate(lab2 = case_when(data_name =="cady_data_ct"~"Any CT definition",
                                         data_name =="cady_data_drug"~"CT drug",
                                         data_name =="cady_data_max_50"~"CT LLN=50 (max)",
                                         data_name =="cady_data_max_53"~"CT LLN=53 (max)",
                                         data_name =="cady_data_mp_50"~"CT LLN=50 (middle point)",
                                         data_name == "cady_data_mp_53"~"CT LLN=53 (middle point)"))
  

  dat <- dat %>% mutate(lab = paste(lab2,"; ",type," (",model,")",sep=""))
  
  dat <- dat %>% filter(marker_name==marker_nam)
  
  dat <- dat %>% arrange(AUC)
  temp <- dat$lab
  dat <- dat %>% mutate(lab = factor(lab,levels = temp))

  
  # Create the forest plot
  if (boncorrect) {
    plot_intervals <- ggplot2::ggplot(dat, ggplot2::aes(x = AUC, y = lab)) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = LB, xmax = UB), 
                              height = 0.2, linewidth = 1.2)+
      # Define labels and titles with large font size
      ggplot2::labs(
        title = paste("Interval Plot of AUC for ",dat$lab1[1],sep=""),
        subtitle = "(With confidence intervals using DeLong method and Bonferroni corrected)",
        x = "Area Under the Curve",
        y = "Setting"
      )
  } else {
    plot_intervals <- ggplot2::ggplot(dat, ggplot2::aes(x = AUC, y = lab)) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = LB_nobon, xmax = UB_nobon), 
                              height = 0.2, linewidth = 1.2)+
      # Define labels and titles with large font size
      ggplot2::labs(
        title = paste("Interval Plot of AUC for ",dat$lab1[1],sep=""),
        subtitle = "(With confidence intervals using DeLong method)",
        x = "Area Under the Curve",
        y = "Setting"
      )
  }
  
  plot_intervals <- plot_intervals +
    ggplot2::geom_point(size = 3) +
    # Add a vertical reference line at 0
    ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed", color = "blue", linewidth = 2)  +
    # Customize the theme to increase font sizes
    ggplot2::theme_minimal(base_size = 18) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 22, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 18),
      axis.title = ggplot2::element_text(size = 20),
      axis.text = ggplot2::element_text(size = 14),
      legend.title = ggplot2::element_text(size = 18),
      legend.text = ggplot2::element_text(size = 16)
    )
  
  plot_intervals
}
