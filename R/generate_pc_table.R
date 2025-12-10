#' @title Generate a Summary Table for a Partly Conditional Cox PH Model
#'
#' @description
#' Creates a publication-ready regression summary table using the \code{gtsummary}
#' package (\code{tbl_regression}) for a partly conditional Cox Proportional Hazards (PH)
#' model (typically fitted with \code{coxph} with a \code{cluster} argument).
#'
#' @param model A fitted regression model object, typically a \code{coxph} object
#'   as provided in the example, or a list containing the model.
#' @param bm_lab A character string specifying the label for the 'marker' variable.
#'   This corresponds to the \code{marker} variable in your example model.
#' @param cov_labs A named list or named vector specifying the labels for the
#'   covariates. The names of the list/vector must match the variable names in
#'   the model, and the values will be the displayed labels.
#' @param title A character string to be used as the caption for the summary table.
#' @param time_lab A character string specifying the label for the 'log_time_to_sample'
#'   variable. Defaults to \code{"log(Time to sample)"}.
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
#' # The function requires a 'coxph' model, typically with a 'cluster' argument
#' # for a "partly conditional" (clustered) Cox model.
#'
#' library(survival)
#' library(gtsummary)
#'
#' # Create a simulated dataset that mimics the structure of your model
#' set.seed(42)
#' N <- 500
#' SubjectID <- rep(1:(N/5), each=5) # 100 subjects
#' data <- data.frame(
#'   t.start = runif(N, 0, 10),
#'   t.stop = runif(N, 10, 50),
#'   status = rbinom(N, 1, 0.1),
#'   log_time_to_sample = runif(N, 0, 5),
#'   marker = rnorm(N, 5, 1.5), # log2(BNP)
#'   Age = runif(N, 40, 80),
#'   lvef_mp_bas = runif(N, 30, 60),
#'   diabetes_mellitus_YN = factor(rbinom(N, 1, 0.3), labels = c("No", "Yes")),
#'   hypertension_YN = factor(rbinom(N, 1, 0.5), labels = c("No", "Yes")),
#'   dyslipidemia_YN = factor(rbinom(N, 1, 0.4), labels = c("No", "Yes")),
#'   treatment_reg = factor(sample(c("Other", "TCH", "Standard"), N, replace = TRUE))
#' )
#' # Correct time variables for Cox model
#' data$t.star <- data$t.stop - data$t.start
#'
#' # Fit the example model (similar to your provided structure)
#' example_fit_pc_cox_BNP <- coxph(
#'   Surv(t.star, status) ~ log_time_to_sample + marker + Age +
#'     lvef_mp_bas + diabetes_mellitus_YN + hypertension_YN +
#'     dyslipidemia_YN + treatment_reg,
#'   data = data,
#'   x = TRUE,
#'   cluster = SubjectID
#' )
#'
#' # --- 2. Define Function Arguments and Call the Function ---
#'
#' # Define the variable labels
#' my_marker_label <- "log2(N-terminal pro B-type natriuretic peptide)"
#' my_covariate_labels <- list(
#'   Age = "Age, years",
#'   lvef_mp_bas = "Baseline LVEF (mid-point), %",
#'   diabetes_mellitus_YN = "Diabetes Mellitus",
#'   hypertension_YN = "Systemic hypertension",
#'   dyslipidemia_YN = "Dyslipidemia",
#'   treatment_reg = "Treatment regime"
#' )
#' my_title <- "Partly Conditional Cox PH Model for Cardiotoxicity (Adjusted by Baseline Characteristics)"
#'
#' # Generate the table
#' table_output <- generate_pc_table(
#'   model = example_fit_pc_cox_BNP,
#'   bm_lab = my_marker_label,
#'   cov_labs = my_covariate_labels,
#'   title = my_title
#' )
#'
#' # Print the resulting table
#' table_output
#' }
#' @export
generate_pc_table <- function(model, bm_lab, cov_labs, title,
                              time_lab = "log(Time to sample)",
                              exponentiate = TRUE, pval_digits = 3) {
  
  # Handle the case where the model is stored inside a list (as in your raw script)
  if (inherits(model, "list") && !inherits(model, "coxph")) {
    stop("The 'model' argument is a list but not a coxph object. Please ensure it is the model object itself or adjust the function to extract it if needed.")
  }
  
  # Combine the fixed labels and the user-provided covariate labels
  # for the 'label' argument of tbl_regression
  all_labels <- c(
    list(log_time_to_sample = time_lab),
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