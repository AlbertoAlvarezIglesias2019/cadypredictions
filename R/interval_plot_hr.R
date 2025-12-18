#' Create a Forest Plot of Hazard Ratios with Confidence Intervals
#'
#' This function generates a forest plot (interval plot) to visualize the
#' Hazard Ratios (HR) and their confidence intervals (CIs) from Cox Proportional
#' Hazards models for a specific biomarker across different study settings or
#' model types.
#'
#' The function relies on the `ggplot2` package for plotting and assumes the
#' input data frame contains specific columns for the HR, confidence interval
#' bounds, and descriptive information about the model and data used.
#'
#' @param dat A data frame containing the results of Cox models. It must include
#'   the following columns:
#'   \itemize{
#'     \item \code{HR}: The estimated Hazard Ratio.
#'     \item \code{LB}: The lower bound of the confidence interval.
#'     \item \code{UB}: The upper bound of the confidence interval.
#'     \item \code{marker_name}: The name of the biomarker (e.g., "BNP", "NT_pro_BNP").
#'     \item \code{data_name}: The name of the dataset/setting used (e.g., "cady_data_ct").
#'     \item \code{type}: The type of analysis (e.g., "Unadjusted", "Adjusted").
#'     \item \code{model}: The specific model name/details.
#'   }
#' @param marker_nam A character string specifying the biomarker to be plotted
#'   (e.g., "BNP", "NT_pro_BNP"). This parameter is used to filter the \code{dat}
#'   and set the main plot title.
#'
#' @return A \code{ggplot} object (a forest plot) visualizing the Hazard Ratios
#'   and their 95\% confidence intervals across different settings.
#'
#' @importFrom dplyr mutate case_when filter arrange
#' @importFrom ggplot2 ggplot aes geom_errorbarh geom_point geom_vline labs theme_minimal theme element_text
#' @export
#'
#' @examples
#' \dontrun{
#' # Create a dummy data frame resembling the expected input structure
#' results_df <- data.frame(
#'   HR = c(1.2, 0.95, 1.5, 0.8),
#'   LB = c(1.05, 0.8, 1.1, 0.65),
#'   UB = c(1.4, 1.1, 2.1, 1.05),
#'   marker_name = c("NT_pro_BNP", "NT_pro_BNP", "NT_pro_BNP", "NT_pro_BNP"),
#'   data_name = c("cady_data_ct", "cady_data_max_50", "cady_data_drug", "cady_data_mp_53"),
#'   type = c("Adjusted", "Unadjusted", "Adjusted", "Unadjusted"),
#'   model = c("Model A", "Model B", "Model C", "Model D")
#' )
#'
#' # Generate the plot for NT-proBNP
#' plot_ntprobnp <- interval_plot_hr(results_df, "NT_pro_BNP")
#' print(plot_ntprobnp)
#' }


interval_plot_hr <- function(dat,marker_nam) {
  
  bc <- 0.05/dim(dat)[1]
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
                                         data_name == "cady_data_mp_53"~"CT LLN=53 (middle point)",
                                         data_name == "cady_data_ct_mace"~"Any CT or MACE (heart failure)",
                                         data_name == "cady_data_ct"~"Any CT",
                                         data_name == "cady_data_mace"~"MACE (heart failure)",
                                         data_name == "cady_data_death"~"Detah from any cause"))
  

  dat <- dat %>% mutate(lab = paste(lab2,"; ",type," (",model,")",sep=""))
  
  dat <- dat %>% filter(marker_name==marker_nam)
  
  dat <- dat %>% arrange(LB)
  temp <- dat$lab
  dat <- dat %>% mutate(lab = factor(lab,levels = temp))

  #####################################
  ## Add bonferroni correction pvalues
  #####################################
  dat <- dat %>% mutate(Significant = if_else(PVAL<bc,"Significant","Not significant"))
  
  # Create the forest plot
  plot_intervals <- ggplot2::ggplot(dat, ggplot2::aes(x = HR, y = lab,color = Significant)) +
    # Add horizontal lines for confidence intervals
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = LB, xmax = UB), 
                            height = 0.2, linewidth = 1.2) +
    # Add points for the estimate
    ggplot2::geom_point(size = 3) +
    # Add a vertical reference line at 0
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "blue", linewidth = 2) +
    # Define labels and titles with large font size
    ggplot2::labs(
      title = paste("Interval Plot of Hazard Ratios for ",dat$lab1[1],sep=""),
      #subtitle = "(Confidence Intervals: Colored by Bonferroni-corrected significance.)",
      x = "Hazard Ratios",
      y = "Setting",
      color = "Bonferroni-adjusted p-values"
    ) +
    # Customize the theme to increase font sizes
    ggplot2::theme_minimal(base_size = 18) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 22, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 18),
      axis.title = ggplot2::element_text(size = 20),
      axis.text = ggplot2::element_text(size = 14),
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 12),
      legend.position = "top"
    )
  
  plot_intervals
}
