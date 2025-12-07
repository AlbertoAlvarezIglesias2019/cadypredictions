#' Create a gtsummary Table Comparing a Single Biomarker by Treatment Regimen
#'
#' This function generates a \code{\link[gtsummary]{tbl_summary}} table for a
#' specific biomarker, comparing its distribution across different treatment
#' regimens. The table displays the **Median [Q1–Q3]** as the primary statistic
#' and, critically, inserts a custom row showing the **count (n)** of unique
#' participants with non-missing data for that biomarker in each treatment group
#' and overall.
#'
#' @param data A data frame in long format. It must include columns named
#'   \code{biomarker}, \code{treatment_reg}, \code{SubjectID}, and the column
#'   specified by \code{value_col}.
#' @param var_name A character string specifying the specific value in the
#'   \code{biomarker} column to summarize (e.g., "CRP", "IL-6").
#' @param value_col A character string specifying the name of the column that
#'   contains the numeric values for the biomarkers (e.g., \code{"Value"}, \code{"Result"}).
#'
#' @return A \code{\link[gtsummary]{tbl_summary}} object that has been styled
#'   with overall column, p-values, and a custom row showing the number of unique
#'   participants with available data for the specified biomarker.
#'
#' @details
#' The function follows a rigorous process to ensure accurate participant counts
#' are displayed alongside the summary statistics:
#'
#' 1.  **Participant Count**: Unique \code{SubjectID}s are counted for the specified
#'     \code{biomarker} and \code{treatment_reg} *only where the value is not missing*.
#' 2.  **Summary Statistic**: The table uses the \code{\link[gtsummary]{tbl_summary}}
#'     function to calculate the **Median [25th–75th Percentile]**.
#' 3.  **Table Modification**: The \code{\link[gtsummary]{add_overall}} and
#'     \code{\link[gtsummary]{add_p}} functions are applied.
#' 4.  **Custom Row Insertion**: A custom row labeled "**Participants (n)**" is
#'     created from the counts calculated in Step 1 and is inserted directly above
#'     the Median statistics row. 
#'
#' This function requires the **gtsummary**, **dplyr**, and **tibble** packages.
#'
#' @importFrom dplyr filter distinct count summarise mutate bind_rows select all_of across
#' @importFrom gtsummary tbl_summary add_overall add_p modify_header all_continuous2 modify_table_body
#' @importFrom tibble tibble add_row
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' # Assuming 'biomarker_data' is a long format dataset
#' # with columns: SubjectID, treatment_reg, biomarker, Value
#'
#' # Summarize 'CRP' biomarker across treatment regimens
#' crp_table <- summarize_biomarker(
#'   data = biomarker_data,
#'   var_name = "CRP",
#'   value_col = "Value"
#' )
#'
#' # Print the gtsummary table
#' print(crp_table)
#' }
#'
summarize_biomarker <- function(data, var_name, value_col) {
  
  # Step 1: Count unique participants per treatment and overall
  counts_tbl <- data %>%
    filter(biomarker == var_name, !is.na(.data[[value_col]])) %>%
    distinct(SubjectID, treatment_reg)
  
  counts_summary <- counts_tbl %>%
    count(treatment_reg, name = "N_participants")
  
  overall_count <- counts_tbl %>%
    summarise(N_participants = dplyr::n_distinct(SubjectID)) %>%
    mutate(treatment_reg = "Overall")
  
  counts_summary <- bind_rows(overall_count, counts_summary) %>%
    mutate(across(N_participants, as.character))
  
  # Step 2: Create summary table
  tbl_main <- data %>%
    filter(biomarker == var_name) %>%
    select(treatment_reg, all_of(value_col)) %>%
    tbl_summary(
      by = treatment_reg,
      type = list(all_of(value_col) ~ "continuous2"),
      statistic = all_continuous2() ~ c(
        "**Median [Q1–Q3]**" = "{median} [{p25} – {p75}]"
      ),
      missing = "no",
      label = list(all_of(value_col) ~ var_name)
    ) %>%
    add_overall() %>%
    add_p() %>%
    modify_header(label = "", all_stat_cols() ~ "**{level}**")
  
  # Step 3: Prepare participant count row
  stat_cols <- grep("^stat_", names(tbl_main$table_body), value = TRUE)
  
  # Initialize empty vector for counts
  counts_row <- as.list(rep(NA_character_, length(stat_cols)))
  names(counts_row) <- stat_cols
  
  # Fill counts in correct order: overall first, then treatment groups
  group_order <- c("Overall", sort(unique(data$treatment_reg)))
  
  for (i in seq_along(group_order)) {
    if (i <= length(stat_cols)) {
      val <- counts_summary$N_participants[counts_summary$treatment_reg == group_order[i]]
      counts_row[[stat_cols[i]]] <- ifelse(length(val) == 1, val, NA_character_)
    }
  }
  
  # Match p.value type
  p_type <- typeof(tbl_main$table_body$p.value)
  p_value_entry <- if (p_type == "double") NA_real_ else NA_character_
  
  # Create counts row tibble
  counts_row <- tibble::tibble(
    variable = var_name,
    row_type = "level",
    label = "Participants (n)",
    !!!counts_row,
    p.value = p_value_entry
  )
  
  # Step 4: Insert counts row safely above the "Median [Q1–Q3]" row
  tbl_main <- tbl_main %>%
    modify_table_body(~{
      body <- .x
      insert_row <- which(grepl("Median", body$label, fixed = TRUE))[1]
      if (is.na(insert_row)) {
        tibble::add_row(body, !!!counts_row, .before = 1)
      } else {
        tibble::add_row(body, !!!counts_row, .before = insert_row)
      }
    })
  
  tbl_main
}