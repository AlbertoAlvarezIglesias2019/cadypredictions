# p-value format for tables with gt_summary

format_pval <- function(x) {
  case_when(
    x < 0.001 ~ "<0.001",
    x > 0.9 ~ ">0.9",
    TRUE ~ sprintf("%.3f", x)  # Format to 3 decimal places
  )
}