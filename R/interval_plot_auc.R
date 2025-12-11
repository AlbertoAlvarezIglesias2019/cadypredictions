#' Create a Forest Plot of Area Under the Curve (AUC) with Confidence Intervals
#'
#' This function generates an interval plot (similar to a forest plot) to visualize the
#' Area Under the Curve (AUC) and its confidence intervals (CIs) for a specific
#' biomarker across different study settings or model types derived from ROC analysis.
#'
#' The function relies on the `ggplot2` package for plotting and assumes the
#' input data frame contains specific columns for the AUC estimate, confidence interval
#' bounds, and descriptive information about the model and data used. The AUC is a
#' common metric in diagnostic testing and indicates the model's ability to
#' discriminate between two classes. 
#'
#' @param dat A data frame containing the results of AUC calculations. It must include
#'   the following columns:
#'   \itemize{
#'     \item \code{AUC}: The estimated Area Under the Curve (AUC).
#'     \item \code{LB}: The lower bound of the confidence interval for the AUC.
#'     \item \code{UB}: The upper bound of the confidence interval for the AUC.
#'     \item \code{marker_name}: The name of the biomarker (e.g., "BNP", "NT_pro_BNP").
#'     \item \code{data_name}: The name of the dataset/setting used (e.g., "cady_data_ct").
#'     \item \code{type}: The type of analysis (e.g., "Unadjusted", "Adjusted").
#'     \item \code{model}: The specific model name/details.
#'   }
#' @param marker_nam A character string specifying the biomarker to be plotted
#'   (e.g., "BNP", "NT_pro_BNP"). This parameter is used to filter the \code{dat}
#'   and set the main plot title. Note the argument name uses \code{marker_nam} as per
#'   the function definition.
#'
#' @return A \code{ggplot} object (an interval plot) visualizing the AUC estimates
#'   and their 95\% confidence intervals across different settings.
#'
#' @importFrom dplyr mutate case_when filter arrange
#' @importFrom ggplot2 ggplot aes geom_errorbarh geom_point geom_vline labs theme_minimal theme element_text
#' @export
#'
#' @examples
#' \dontrun{
#' # Create a dummy data frame resembling the expected input structure
#' results_auc_df <- data.frame(
#'   AUC = c(0.75, 0.82, 0.68, 0.90),
#'   LB = c(0.68, 0.75, 0.60, 0.85),
#'   UB = c(0.82, 0.89, 0.76, 0.95),
#'   marker_name = c("NT_pro_BNP", "NT_pro_BNP", "NT_pro_BNP", "NT_pro_BNP"),
#'   data_name = c("cady_data_ct", "cady_data_max_50", "cady_data_drug", "cady_data_mp_53"),
#'   type = c("Adjusted", "Unadjusted", "Adjusted", "Unadjusted"),
#'   model = c("Model A", "Model B", "Model C", "Model D")
#' )
#'
#' # Generate the plot for NT-proBNP
#' plot_auc_ntprobnp <- interval_plot_auc(results_auc_df, "NT_pro_BNP")
#' print(plot_auc_ntprobnp)
#' }
#' 
#' 


interval_plot_auc <- function(dat,marker_nam) {
  
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
  plot_intervals <- ggplot2::ggplot(dat, ggplot2::aes(x = AUC, y = lab)) +
    # Add horizontal lines for confidence intervals
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = LB, xmax = UB), 
                            height = 0.2, linewidth = 1.2) +
    # Add points for the estimate
    ggplot2::geom_point(size = 3) +
    # Add a vertical reference line at 0
    ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed", color = "blue", linewidth = 2) +
    # Define labels and titles with large font size
    ggplot2::labs(
      title = paste("Interval Plot of AUC for ",dat$lab1[1],sep=""),
      subtitle = "(With confidence intervals using DeLong method)",
      x = "Area Under the Curve",
      y = "Setting"
    ) +
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
