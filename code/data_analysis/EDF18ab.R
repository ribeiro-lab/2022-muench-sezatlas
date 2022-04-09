# plot -----------------------------------------------------------------


data <- read.csv('57C10_flyPAD.csv')

plot_c <-
  function(i,
           p_y_y,
           p_s_y,
           ytransform = waiver(),
           ylab,
           y.lim = NA,
           cex = 5) {
    plot_data <-
      filter(data, condition %in% c("mated_10dS", "mated_ff"), id == i) %>%
      mutate(
        condition = fct_recode(
          condition,
          "fully fed" = "mated_ff",
          "protein\ndeprived" = "mated_10dS"
        ),
        condition = fct_relevel(condition, "fully fed")
      )
    
    
    stat_wtest <- plot_data %>%
      group_by(spot, experiment) %>%
      wilcox_test(value ~ condition) %>%
      add_significance()
    
    
    p <-
      ggplot(data = filter(plot_data, spot == 'yeast'), aes(x = condition, y = value)) +
      ggbeeswarm::geom_beeswarm(
        aes(color = condition),
        size = 3.5,
        shape = 16,
        cex = cex
      ) +
      geom_boxplot(
        color = 'grey40',
        alpha = 0,
        width = .5,
        outlier.alpha = 0,
        size = 1,
        outlier.shape = 16,
        outlier.size = 3.5,
        outlier.colour = "black"
      ) +
      scale_color_manual(values = fs_colors2) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)), labels =
                           ytransform) +
      theme_dahaniel2(base_size = 24) +
      labs(x = NULL, y =  ylab[1], title = 'yeast') +
      coord_cartesian(clip = 'off') +
      stat_pvalue_manual(
        filter(stat_wtest, spot == 'yeast'),
        label = "{formatC(p,2)}",
        size = 5,
        bracket.size = .6,
        tip.length = 0,
        y.position = c(p_y_y)
      ) +
      theme(legend.position = 'none',
            panel.border = element_blank()) +
      expand_limits(y = 0)
    if (!is.na(y.lim[1]))
      p  <- p + coord_cartesian(ylim = y.lim[c(1, 2)])
    ggsave(
      paste0(out_path, 'Fig6_bp_wt_', i, '_yeast.pdf'),
      width = 25,
      height = 45,
      dpi = 300,
      scale = 3,
      units = 'mm',
      plot = p
    )
    
    p <-
      ggplot(data = filter(plot_data, spot == 'sucrose'),
             aes(x = condition, y = value)) +
      ggbeeswarm::geom_beeswarm(
        aes(color = condition),
        size = 3.5,
        shape = 16,
        cex = cex
      ) +
      geom_boxplot(
        color = 'grey40',
        alpha = 0,
        width = .5,
        outlier.alpha = 0,
        size = 1,
        outlier.shape = 16,
        outlier.size = 3.5,
        outlier.colour = "black"
      ) +
      scale_color_manual(values = fs_colors2) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)), labels =
                           ytransform) +
      theme_dahaniel2(base_size = 24) +
      labs(x = NULL, y =  ylab[2], title = 'sucrose') +
      coord_cartesian(clip = 'off') +
      stat_pvalue_manual(
        filter(stat_wtest, spot == 'sucrose'),
        label = "{formatC(p,2)}",
        size = 5,
        bracket.size = .6,
        tip.length = 0,
        y.position = c(p_s_y)
      ) +
      theme(legend.position = 'none',
            panel.border = element_blank()) +
      expand_limits(y = 0)
    if (!is.na(y.lim[1]))
      p  <- p + coord_cartesian(ylim = y.lim[c(3, 4)])
    ggsave(
      paste0(out_path, 'Fig6_bp_wt_', i, '_sucrose.pdf'),
      width = 25,
      height = 45,
      dpi = 300,
      scale = 3,
      units = 'mm',
      plot = p
    )
  }


plot_c('FB_IBI', 120,400, y.lim = c(0,120,0,400), ylab = c('mean inter burst interval [s]', 'mean inter burst interval [s]'), cex = 2.5)
plot_c('sips_per_burst', 22,20, ylab = c('sips per burst', 'sips per burst'))








