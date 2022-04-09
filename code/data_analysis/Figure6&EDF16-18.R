# plot borboleta overlap ----------------------------------------------
ovrlp <- read.csv('braincode_overlap.csv')
ovrlp <- ovrlp %>% filter(region_bc != 0) %>%
  mutate(overlap_mp = overlap / size_mp,
         overlap_bc = overlap / size_bc)

ggplot(filter(ovrlp, region_mp == "(21, 9)")) +
  geom_bar(aes(x = reorder(region_bc, overlap), y = overlap), stat = 'identity', fill = '#d96fc9ff') +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = 'braincode cluster', y = 'overlapping voxels', title = 'overlap with outer borboleta') +
  theme_dahaniel2() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_blank()
  )
ggsave(
  paste0(out_path, "bp_braincode_overlap_outer-borboleta.pdf"),
  width = 30,
  height = 20,
  dpi = plot_dpi,
  scale = plot_scale,
  units = plot_units,
  plot = last_plot()
)
write.csv(
  filter(ovrlp, region_mp == "(21, 9)") %>% select(region_mp, region_bc, overlap),
  paste0(out_path, "source-data_Fig6b_outer.csv")
)

ggplot(filter(ovrlp, region_mp == "(57, 60)")) +
  geom_bar(aes(x = reorder(region_bc, overlap), y = overlap), stat = 'identity', fill = '#b5d96fff') +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = 'braincode cluster', y = 'overlapping voxels', title = 'overlap with inner borboleta') +
  theme_dahaniel2() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_blank()
  )
ggsave(
  paste0(out_path, "bp_braincode_overlap_inner-borboleta.pdf"),
  width = 30,
  height = 20,
  dpi = plot_dpi,
  scale = plot_scale,
  units = plot_units,
  plot = last_plot()
)
write.csv(
  filter(ovrlp, region_mp == "(57, 60)") %>% select(region_mp, region_bc, overlap),
  paste0(out_path, "source-data_Fig6b_inner.csv")
)




# plot flyPAD data --------------------------------------------------------
data <- read.csv("TrpA1_flyPAD.csv")
# label categories
data$borboleta <- 'inner'
data$borboleta[which(data$condition %in% c('22B10', '78E06', '78G03'))] <-
  'outer'
data$borboleta[which(data$condition %in% c('empty'))] <- 'control'
borbo_col <-
  c('inner' = '#B5D96F',
    'outer' = '#D96FC9',
    'control' = 'grey70')

formatter1000 <- function(x)
  x / 1000

plot_c <-
  function(i,
           c,
           p_y_y,
           p_s_y,
           ytransform = waiver(),
           ylab,
           y.lim = NA,
           cex = 5) {
    c_exp <- unique(filter(data, condition == c)$experiment)
    plot_data <-
      filter(data, experiment %in% c_exp, condition %in% c(c, 'empty'), id == i)
    
    stat_wtest <- plot_data %>%
      group_by(spot, experiment) %>%
      wilcox_test(value ~ condition, ref.group = 'empty') %>%
      add_significance()
    
    plot_data$condition <-
      factor(plot_data$condition, levels = sort(unique(plot_data$condition), decreasing = TRUE))
    p <-
      ggplot(data = filter(plot_data, spot == 'yeast'), aes(x = condition, y = value)) +
      ggbeeswarm::geom_beeswarm(
        aes(color = borboleta),
        size = .5,
        shape = 16,
        cex = cex
      ) +
      geom_boxplot(
        color = 'grey40',
        alpha = 0,
        width = .5,
        size = .3,
        outlier.shape = NA
      ) +
      scale_color_manual(values = borbo_col) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.07)), labels =
                           ytransform) +
      theme_dahaniel2() +
      labs(x = NULL, y =  ylab[1]) +
      coord_cartesian(clip = 'off') +
      stat_pvalue_manual(
        filter(stat_wtest, spot == 'yeast'),
        label = "{formatC(p, 2)}",
        size = stat_size,
        bracket.size = stat_bracket.size,
        tip.length = stat_tip.length,
        y.position = c(p_y_y)
      ) +
      theme(
        legend.position = 'none',
        panel.border = element_blank(),
        axis.ticks.x = element_blank()
      ) +
      expand_limits(y = 0)
    if (!is.na(y.lim[1]))
      p  <- p + coord_cartesian(ylim = y.lim[c(1, 2)])
    ggsave(
      paste0(out_path, 'Fig6_bp_behavior_', c, '_', i, '_yeast.pdf'),
      width = 15,
      height = 28,
      dpi = plot_dpi,
      scale = plot_scale,
      units = plot_units,
      plot = p
    )
    write.csv(
      filter(plot_data, spot == 'yeast') %>% select(spot, id, condition, value),
      paste0(out_path, "/source-data_behavior_", c, "_", i, "_yeast.csv")
    )
    print(
      paste(
        'median yeast:',
        paste(layer_data(p, 2)$middle, collapse = ", "),
        "delta: ",
        layer_data(p, 2)$middle[2] - layer_data(p, 2)$middle[1]
      )
    )
    print(layer_data(p, 2))
    
    p <-
      ggplot(data = filter(plot_data, spot == 'sucrose'),
             aes(x = condition, y = value)) +
      ggbeeswarm::geom_beeswarm(
        aes(color = borboleta),
        size = .5,
        shape = 16,
        cex = cex
      ) +
      geom_boxplot(
        color = 'grey40',
        alpha = 0,
        width = .5,
        size = .3,
        outlier.shape = NA
      ) +
      scale_color_manual(values = borbo_col) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.07)), labels =
                           ytransform) +
      theme_dahaniel2() +
      labs(x = NULL, y =  ylab[2]) +
      coord_cartesian(clip = 'off') +
      stat_pvalue_manual(
        filter(stat_wtest, spot == 'sucrose'),
        label = "{formatC(p, 2)}",
        size = stat_size,
        bracket.size = stat_bracket.size,
        tip.length = stat_tip.length,
        y.position = c(p_s_y)
      ) +
      theme(
        legend.position = 'none',
        panel.border = element_blank(),
        axis.ticks.x = element_blank()
      ) +
      expand_limits(y = 0)
    if (!is.na(y.lim[1]))
      p  <- p + coord_cartesian(ylim = y.lim[c(3, 4)])
    ggsave(
      paste0(out_path, 'Fig6_bp_behavior_', c, '_', i, '_sucrose.pdf'),
      width = 15,
      height = 28,
      dpi = plot_dpi,
      scale = plot_scale,
      units = plot_units,
      plot = p
    )
    write.csv(
      filter(plot_data, spot == 'sucrose') %>% select(spot, id, condition, value),
      paste0(out_path, "/source-data_behavior_", c, "_", i, "_sucrose.csv")
    )
    print(paste('median sucrose:', paste(layer_data(p, 2)$middle, collapse = ", ")))
  }

plot_c('sips_n', '78E06', 4800,1250, y.lim = c(0,4800,0,1250), ytransform = formatter1000, ylab = c('sips on yeast (×1000)', 'sips on sucrose (×1000)'))
plot_c('sips_n', '78G03', 4800,1250, y.lim = c(0,4800,0,1250), ytransform = formatter1000, ylab = c('sips on yeast (×1000)', 'sips on sucrose (×1000)'))
plot_c('sips_n', '22B10', 4800,1250, y.lim = c(0,4800,0,1250), ytransform = formatter1000, ylab = c('sips on yeast (×1000)', 'sips on sucrose (×1000)'))
plot_c('sips_n', '49E02', 4800,1250, y.lim = c(0,4800,0,1250), ytransform = formatter1000, ylab = c('sips on yeast (×1000)', 'sips on sucrose (×1000)'))
plot_c('sips_n', '70C11', 4800,1250, y.lim = c(0,4800,0,1250), ytransform = formatter1000, ylab = c('sips on yeast (×1000)', 'sips on sucrose (×1000)'))
plot_c('sips_n', '77C10', 4800,1250, y.lim = c(0,4800,0,1250), ytransform = formatter1000, ylab = c('sips on yeast (×1000)', 'sips on sucrose (×1000)'))

plot_c('sips_per_burst', '78E06', 24,12, ylab = c('sips per burst', 'sips per burst'))
plot_c('sips_per_burst', '78G03', 18,12.8, ylab = c('sips per burst', 'sips per burst'))
plot_c('sips_per_burst', '22B10', 25,12, ylab = c('sips per burst', 'sips per burst'))
plot_c('sips_per_burst', '49E02', 19,13, ylab = c('sips per burst', 'sips per burst'))
plot_c('sips_per_burst', '70C11', 15.5,16.2, ylab = c('sips per burst', 'sips per burst'))
plot_c('sips_per_burst', '24E08', 24.5,12, ylab = c('sips per burst', 'sips per burst'))
plot_c('sips_per_burst', '77C10', 13.5,13.2, ylab = c('sips per burst', 'sips per burst'))

plot_c('FB_IBI', '78E06', 200,200, ylab = c('mean inter burst interval [s]', 'mean inter burst interval [s]'), y.lim = c(0,200,0,200), cex = 3.5)
plot_c('FB_IBI', '78G03', 160,420, ylab = c('mean inter burst interval [s]', 'mean inter burst interval [s]'), y.lim = c(0,170,0,420), cex = 3.5)
plot_c('FB_IBI', '22B10', 185,165, ylab = c('mean inter burst interval [s]', 'mean inter burst interval [s]'), y.lim = c(0,185,0,180), cex = 3.5)
plot_c('FB_IBI', '49E02', 250,250, ylab = c('mean inter burst interval [s]', 'mean inter burst interval [s]'), y.lim = c(0,250,0,250), cex = 3.5)
plot_c('FB_IBI', '70C11', 255,240, ylab = c('mean inter burst interval [s]', 'mean inter burst interval [s]'), y.lim = c(0,255,0,250), cex = 3.5)
plot_c('FB_IBI', '24E08', 200,300, ylab = c('mean inter burst interval [s]', 'mean inter burst interval [s]'), y.lim = c(0,200,0,300), cex = 3.5)
plot_c('FB_IBI', '77C10', 150,340, ylab = c('mean inter burst interval [s]', 'mean inter burst interval [s]'), y.lim = c(0,150,0,340), cex = 3.5)






# tracking micromovements -----------------------------------------------------
data_tracking <- read.csv('TrpA1_tracking_yeast_micromovements.csv')

data_tracking <- data_tracking %>%
  mutate(
    condition = recode(condition, Empty = 'empty'),
    condition = fct_relevel(condition, 'empty', after = 0)
  ) %>%
  mutate(
    borboleta = 'outer',
    borboleta = if_else(condition == 'empty', 'control', 'outer')
  ) %>%
  group_by(fly, condition, borboleta) %>%
  summarize(
    visits = n(),
    total_duration = sum(duration) / 60,
    mean_duration = mean(duration)
  ) %>%
  group_by(condition) %>%
  mutate(n = n(), spot = 'yeast') %>%
  ungroup()

stat_wtest <- data_tracking %>%
  group_by(spot) %>%
  wilcox_test(total_duration ~ condition, ref.group = 'empty') %>%
  add_significance()

ggplot(data = filter(data_tracking, spot == 'yeast'),
       aes(x = condition, y = total_duration)) +
  ggbeeswarm::geom_beeswarm(
    aes(color = borboleta),
    size = 2.8,
    shape = 16,
    alpha = 1,
    cex = 2
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
  scale_color_manual(values = borbo_col) +
  theme_dahaniel2(base_size = 24) +
  scale_y_continuous(expand = expansion(c(0.05, 0.07)), breaks = c(0, 5, 10)) +
  labs(x = 'Gal4 line', y =  'total duration of\nmicromovements on yeast [min]') +
  stat_pvalue_manual(
    stat_wtest,
    label = "{formatC(p,2)}",
    size = 5,
    bracket.size = .6,
    tip.length = 0.01,
    y.position = c(10.8, 11.8)
  ) +
  coord_cartesian(clip = 'on', ylim = c(0, 11.5)) +
  theme(legend.position = 'none',
        panel.border = element_blank())
ggsave(
  paste0(out_path, "bp_tracking-mm_Fig6_yeast.svg"),
  width = 35,
  height = 60,
  dpi = 300,
  scale = 3,
  units = 'mm',
  plot = last_plot()
)

write.csv(
  filter(data_tracking, spot == 'yeast') %>% select(spot, condition, total_duration),
  paste0(out_path, "source-data_bp_tracking-mm_Fig6_yeast.csv")
)