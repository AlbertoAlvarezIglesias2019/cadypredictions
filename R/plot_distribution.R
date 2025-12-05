
plot_distribution <- function(data, lab_variable, visit, title, binwidth) {
  # Filter for the specific visit and lab variable
  disp_summary <- data %>%
    filter(Visit == visit, Variable == lab_variable) %>%
    group_by(RandArm) %>%
    summarize(mean_val = mean(Value, na.rm = TRUE),
              median_val = median(Value, na.rm = TRUE), .groups = "drop")
  
  # Generate the histogram and density plot
  plot <- data %>%
    filter(Visit == visit, Variable == lab_variable) %>%
    ggplot(aes(x = Value)) +
    geom_histogram(aes(y = ..density..), binwidth = binwidth, fill = "skyblue", color = "black", alpha = 0.5) +
    geom_density(color = "black", size = 0.5) +
    geom_vline(data = disp_summary, aes(xintercept = mean_val, color = "Mean"), linetype = "dashed", size = 1) +
    geom_vline(data = disp_summary, aes(xintercept = median_val, color = "Median"), linetype = "dashed", size = 1) +
    facet_wrap(~ RandArm) +
    theme_minimal() +
    labs(title = title, x = lab_variable, y = "Density") +
    scale_color_manual(values = c("Mean" = "blue", "Median" = "red"))
  
  return(plot)
}