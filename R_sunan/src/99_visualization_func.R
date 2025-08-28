plot_pw_effects <- function(dataframe,
                            title = "Point-wise true effects throughout the day",
                            xlab = "Time of Day (Hours)",
                            ylab = expression(beta),
                            facet_labels = c(
                              "Beta" = "beta[Age]",
                              "Intercept" = "beta[0]"),
                            type_vars = c("Beta", 
                                          "Intercept")) {
  
  dataframe |>
    mutate(Minutes = as.numeric(sub('min_', '', Variable)),
           Hour = Minutes / 60) |>
    dplyr::select(-Variable, -Minutes) |>
    tidyr::pivot_longer(
      cols = all_of(type_vars),
      names_to = "Type",
      values_to = "Value"
    ) |>
    ggplot(aes(x = Hour, y = Value)) +
    geom_point(color = "blue", alpha = 0.6, size = 0.8) +
    facet_wrap(~ Type,
               scales = "free_y",
               labeller = as_labeller(facet_labels, default = label_parsed)) +
    scale_x_continuous(
      breaks = seq(0, 24, by = 6),
      labels = function(x) paste0(x, ":00")
    ) +
    labs(title = title, x = xlab, y = ylab) +
    theme(
      text = element_text(size = 14),
      axis.title.x = element_text(face = 'bold'),
      axis.title.y = element_text(face = 'bold'),
      strip.text = element_text(face = "bold", size = 12)
    )
}


visualize_from_folder <- function(folder = NA,
                                  file_pattern = '',
                                  L = NA){
  all_results <- list.files(path = folder, 
                            pattern = file_pattern, 
                            full.names = TRUE) %>% map_dfr(readRDS)
  
  summary_design = all_results %>%
    group_by(boot_type, var) %>%
    summarise(
      MISE = mean(MISE),
      mean_joint_se = mean(mean_joint_se),
      mean_pw_se = mean(mean_pw_se),
      cover_joint_ratio = mean(cover_joint),
      cover_pw_ratio = mean(cover_pw) / L,
      .groups = 'drop'
    )
  
  formatted_table <- summary_design %>%
    mutate(across(where(is.numeric), ~round(., 4)))  
  kable(formatted_table, digits = 4,
        caption = "Summary Statistics", format = "html")
  # return(formatted_table)
}
