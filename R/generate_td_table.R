#' @title Generate a Summary Table for a Time-Dependent Covariates Cox PH Model
#'
#' @description
#' Creates a publication-ready regression summary table using the \code{gtsummary}
#' package (\code{tbl_regression}) for a Cox Proportional Hazards (PH) model
#' that includes time-dependent covariates (TDC).
#'
#' @param model A fitted regression model object, typically a \code{coxph} object
#'   from a time-dependent analysis.
#' @param bm_lab A character string specifying the label for the 'marker' variable
#'   (the time-dependent covariate in this context).
#' @param cov_labs A named list or named vector specifying the labels for the
#'   baseline covariates. The names of the list/vector must match the variable
#'   names in the model, and the values will be the displayed labels.
#' @param title A character string to be used as the caption for the summary table.
#' @param exponentiate A logical indicating whether to exponentiate the coefficient estimates
#'   (e.g., convert log-Hazard Ratios to Hazard Ratios). Defaults to \code{TRUE}.
#' @param pval_digits An integer specifying the number of decimal places for p-values.
#'   Defaults to \code{3}.
#'
#' @return A \code{gtsummary} object (\code{tbl_regression} class) that can be
#'   printed directly or further modified.
#'
#' @import gtsummary
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' # --- 1. Prepare Dummy Data and Model for Example ---
#' # The function requires a 'coxph' model, which can be fit using data
#' # prepared with the 'tmerge' function for a TDC analysis.
#'
#' library(survival)
#' library(gtsummary)
#'
#' # Create a simulated dataset that mimics the structure of a TDC analysis
#' set.seed(42)
#' N <- 500
#' data_long <- data.frame(
#'   id = 1:N,
#'   tstart = 0,
#'   tstop = rexp(N, 1/30),
#'   status = rbinom(N, 1, 0.1),
#'   marker = rnorm(N, 5, 1.5), # This will be treated as the TDC
#'   Age = runif(N, 40, 80),
#'   lvef_mp_bas = runif(N, 30, 60),
#'   diabetes_mellitus_YN = factor(rbinom(N, 1, 0.3), labels = c("No", "Yes")),
#'   hypertension_YN = factor(rbinom(N, 1, 0.5), labels = c("No", "Yes")),
#'   dyslipidemia_YN = factor(rbinom(N, 1, 0.4), labels = c("No", "Yes")),
#'   treatment_reg = factor(sample(c("Other", "TCH", "Standard"), N, replace = TRUE))
#' )
#'
#' # Fit the example Time-Dependent Covariates model
#' # (A simple example without true tmerge, just to represent the formula structure)
#' example_fit_td_cox_BNP <- coxph(
#'   Surv(tstart, tstop, status) ~ marker + Age + lvef_mp_bas +
#'     diabetes_mellitus_YN + hypertension_YN + dyslipidemia_YN +
#'     treatment_reg,
#'   data = data_long
#' )
#'
#' # --- 2. Define Function Arguments and Call the Function ---
#'
#' # Define the variable labels
#' my_marker_label <- "log2(BNP) (Time-dependent)"
#' my_covariate_labels <- list(
#'   Age = "Age, years",
#'   lvef_mp_bas = "Baseline LVEF (mid-point), %",
#'   diabetes_mellitus_YN = "Diabetes Mellitus",
#'   hypertension_YN = "Systemic hypertension",
#'   dyslipidemia_YN = "Dyslipidemia",
#'   treatment_reg = "Treatment regime"
#' )
#' my_title <- "Time-Dependent Covariates Cox PH Model for Cardiotoxicity"
#'
#' # Generate the table
#' table_output <- generate_td_table(
#'   model = example_fit_td_cox_BNP,
#'   bm_lab = my_marker_label,
#'   cov_labs = my_covariate_labels,
#'   title = my_title
#' )
#'
#' # Print the resulting table
#' table_output
#' }
#' @export
generate_td_table <- function(model, bm_lab, cov_labs, title,
                              exponentiate = TRUE, pval_digits = 3) {
  
  # Combine the marker label and the user-provided covariate labels
  # for the 'label' argument of tbl_regression
  all_labels <- c(
    list(marker = bm_lab),
    cov_labs
  )
  
  # 1. Create the regression table using tbl_regression
  tbl_output <- tbl_regression(
    x = model,
    exponentiate = exponentiate,
    pvalue_fun = label_style_pvalue(digits = pval_digits),
    label = all_labels
  )
  
  # 2. Modify the header to remove the 'Characteristic' or 'Variable' label
  tbl_output <- tbl_output %>%
    modify_header(label ~ "")
  
  # 3. Add the user-defined caption (title)
  tbl_output <- tbl_output %>%
    modify_caption(title)
  
  # Return the gtsummary object
  return(tbl_output)
}