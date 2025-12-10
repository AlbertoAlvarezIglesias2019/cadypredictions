#' @title Fit Cox PH Model and Generate gtsummary Regression Table
#'
#' @description This function fits a Cox Proportional Hazards model and
#'              generates a clean, formatted regression table using gtsummary,
#'              allowing for easy swapping of the biomarker and its label.
#'
#' @param data A data frame containing the variables for the Cox model.
#'             The data should already be filtered/subsetted as needed (e.g., to the 'training' set).
#' @param time_event_vars A character vector of length 2 specifying the column names for
#'                        the survival time and event status, respectively.
#'                        (e.g., c("time_to_event", "status")).
#' @param covariates_coxph A character vector of covariates (column names) to include
#'                         in the model, excluding the biomarker.
#' @param biomarker A character string representing the right-hand side (RHS) expression
#'                  for the biomarker, often including a transformation (e.g., "log2(BNP_bl)").
#' @param bm_lab A character string: the label for the selected biomarker (e.g., "log2(NT-proBNP)").
#' @param cov_labs A named list of labels for all non-biomarker covariates (e.g., Age, LVEF, etc.).
#'                 The names must match the variable names in `covariates_coxph`.
#' @param title A character string to be used as the caption for the gtsummary table.
#' @param ... Additional arguments passed to the `coxph` function (e.g., `ties`).
#'
#' @return A list containing two elements:
#'         \itemize{
#'           \item \code{coxph_model}: The fitted \code{coxph} model object.
#'           \item \code{tbl_regression}: The \code{tbl_regression} object from gtsummary.
#'         }
#'
#' @export
generate_coxph_table <- function(data, time_event_vars = c("time_to_event", "status"),
                                 covariates_coxph, biomarker, bm_lab, cov_labs, title, ...) {
  
  # 1. Input Validation (basic check)
  if (!requireNamespace("survival", quietly = TRUE) || !requireNamespace("gtsummary", quietly = TRUE)) {
    stop("Packages 'survival' and 'gtsummary' are required. Please install them.")
  }
  
  # 2. Prepare the combined labels list
  # Create the biomarker label list item dynamically
  biomarker_lab_list <- list(bm_lab)
  names(biomarker_lab_list) <- biomarker
  
  # Combine the biomarker label with the static covariate labels
  all_model_labels <- c(biomarker_lab_list, cov_labs)
  
  # 3. Build the right-hand side (RHS) of the model
  rhs <- c(biomarker, covariates_coxph)
  
  # 4. Build the formula: Surv(time, status) ~ biomarker + cov1 + cov2 + ...
  response_formula <- paste0("Surv(", time_event_vars[1], ", ", time_event_vars[2], ")")
  fml <- reformulate(rhs, response = response_formula)
  
  # 5. Fit the Cox PH model
  fit_coxph <- coxph(fml, data = data, ...)
  
  # 6. Clean and format the table summary
  tbl_fit_coxph <- tbl_regression(
    fit_coxph,
    exponentiate = TRUE,
    pvalue_fun = label_style_pvalue(digits = 3),
    label = all_model_labels # Use the combined, correct labels list
  ) %>%
    modify_header(label ~ "") %>%
    modify_caption(title)
  
  # 7. Return the desired outputs
  return(
    list(
      coxph_model = fit_coxph,
      tbl_regression = tbl_fit_coxph
    )
  )
}