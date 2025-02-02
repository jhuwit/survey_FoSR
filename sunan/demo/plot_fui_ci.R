plot_fui_ci <- function(dat,
                        cross_interval = TRUE,
                        label = 'please enter your figure label',
                        sublabel = NA,
                        ylim = NA,
                        yinterval = NA){
  if (is.na(ylim[1])) {
    ylim = c(min(dat$CI.lower.joint), max(dat$CI.upper.joint))
    yinterval = round((max(dat$CI.upper.joint)-min(dat$CI.upper.joint))/6)
  }

  colnames(dat) = c('time', 'beta',
                    'lower.pointwise', 'upper.pointwise',
                    'lower.joint', 'upper.joint')

  dat$joint.sig = ifelse(sign(dat$lower.joint) == sign(dat$upper.joint), 1,0)
  dat$point.sig = ifelse(sign(dat$lower.pointwise) == sign(dat$upper.pointwise), 1,0)
  dat$label.sig = ifelse(dat$joint.sig == 1, 2, dat$point.sig)

  if (cross_interval == TRUE) {
    ## Build Across Confidence Interval
    cross_intervals <- dat %>%
      mutate(pointwise_cross = ifelse((lower.pointwise < 0 & upper.pointwise < 0) |
                                        (lower.pointwise > 0 & upper.pointwise > 0), TRUE, FALSE),
             joint_cross = ifelse((lower.joint < 0 & upper.joint < 0) |
                                    (lower.joint > 0 & upper.joint > 0), TRUE, FALSE)) %>%
      mutate(group = cumsum(pointwise_cross != lag(pointwise_cross,
                                                   default = first(pointwise_cross)) |
                              joint_cross != lag(joint_cross,
                                                 default = first(joint_cross)))) %>%
      group_by(group) %>%
      summarise(start = min(time) / 60,
                end = max(time) / 60,
                pointwise_cross = first(pointwise_cross),
                joint_cross = first(joint_cross)) %>%
      filter(pointwise_cross | joint_cross)


    plt = ggplot(dat, aes(x = time / 60, y = beta)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      geom_ribbon(aes(ymin = lower.joint, ymax = upper.joint),
                  fill = "#d6d6d6", alpha = 0.5) +
      geom_ribbon(aes(ymin = lower.pointwise, ymax = upper.pointwise),
                  fill = "#8b8b8b", alpha = 0.5) +

      geom_rect(data = cross_intervals, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf,
                                            fill = ifelse(joint_cross, "Joint CI", "Pointwise CI")),
                alpha = 0.1, inherit.aes = FALSE) +
      scale_fill_manual(values = c("Pointwise CI" = "red", "Joint CI" = "blue")) +
      labs(title = label, x = "Time (hour)", y = 'Physical Intensity (count/min)',
           fill = "Significance") +
      geom_line(aes(y = beta), color = "black") +
      scale_x_continuous(breaks = c(1,6,12,18,23),
                         labels = function(x) paste0(sprintf("%02d", x), ":00")) +
      scale_y_continuous(breaks = seq(ylim[1], ylim[2], by = yinterval), limits = ylim) +
      theme(
        panel.background = element_blank(),
        # axis.ticks = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        axis.line = element_line(size = 0.5),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
      ) +
      theme(legend.text = element_text(size = 12, face = 'bold'),
            legend.title = element_text(size = 12, face = 'bold'),
            plot.title = element_text(size = 12),
            axis.text = element_text(size = 10),
            axis.title = element_text(size = 10, face = "bold"),
            panel.grid.major = element_line(color = "gray90", size = 0.2),
            panel.grid.minor = element_blank()) +
      guides(color = guide_legend(title = NULL), linetype = guide_legend(title = NULL))
  }
  if (cross_interval == FALSE){
    plt = ggplot(dat, aes(x = time / 60, y = beta)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      geom_ribbon(aes(ymin = lower.joint, ymax = upper.joint),
                  fill = "#d6d6d6", alpha = 0.5) +
      geom_ribbon(aes(ymin = lower.pointwise, ymax = upper.pointwise),
                  fill = "#8b8b8b", alpha = 0.5) +
      scale_fill_manual(values = c("Pointwise CI" = "red", "Joint CI" = "blue")) +
      labs(title = label,
           subtitle = sublabel,
           x = "Time (hour)", y = 'Physical Intensity (count/min)',
           fill = "Significance") +
      geom_line(aes(y = beta), color = "black") +
      scale_x_continuous(breaks = c(1,6,12,18,23),
                         labels = function(x) paste0(sprintf("%02d", x), ":00")) +
      scale_y_continuous(breaks = seq(ylim[1], ylim[2], by = yinterval), limits = ylim) +
      theme(
        panel.background = element_blank(),
        # axis.ticks = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        axis.line = element_line(size = 0.5),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
      ) +
      theme(legend.text = element_text(size = 12, face = 'bold'),
            legend.title = element_text(size = 12, face = 'bold'),
            plot.title = element_text(size = 12),
            axis.text = element_text(size = 10),
            axis.title = element_text(size = 10, face = "bold"),
            panel.grid.major = element_line(color = "gray90", size = 0.2),
            panel.grid.minor = element_blank()) +
      guides(color = guide_legend(title = NULL), linetype = guide_legend(title = NULL))
  }

  return(plt)
}
