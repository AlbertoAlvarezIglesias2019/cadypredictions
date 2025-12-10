#' @title Fit Cox PH Model and Generate gtsummary Regression Table
#'
#' @description This function fits a Cox Proportional Hazards model and
#'              generates a clean, formatted regression table using gtsummary.
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
#' @param labels A named list or named vector where names are the variable names/expressions
#'               in the model and values are the desired labels for the final table.
#'               (e.g., list("log2(BNP_bl)" = "log2(BNP)", "Age" = "Age, years")).
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
generate_coxph <- function(data, time_event_vars = c("time_to_event", "status"),
                                 covariates_coxph, biomarker, labels, title, ...) {
  
  # 1. Input Validation (basic check for required packages)
  if (!requireNamespace("survival", quietly = TRUE) || !requireNamespace("gtsummary", quietly = TRUE)) {
    stop("Packages 'survival' and 'gtsummary' are required. Please install them.")
  }
  
  # 2. Build the right-hand side (RHS) of the model
  # The biomarker expression goes first, followed by the standard covariates
  rhs <- c(biomarker, covariates_coxph)
  
  # 3. Build the formula: Surv(time, status) ~ biomarker + cov1 + cov2 + ...
  response_formula <- paste0("Surv(", time_event_vars[1], ", ", time_event_vars[2], ")")
  fml <- reformulate(rhs, response = response_formula)
  
  # 4. Fit the Cox PH model
  fit_coxph <- coxph(fml, data = data, ...)
  
  # 5. Clean and format the table summary
  tbl_fit_coxph <- tbl_regression(
    fit_coxph,
    exponentiate = TRUE,
    pvalue_fun = label_style_pvalue(digits = 3),
    label = labels # Use the provided labels list
  ) %>%
    modify_header(label ~ "") %>%
    modify_caption(title)
  
  # 6. Return the desired outputs in a named list
  return(
    list(
      coxph_model = fit_coxph,
      tbl_regression = tbl_fit_coxph
    )
  )
}