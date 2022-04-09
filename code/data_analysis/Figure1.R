data_fp <- read.csv('57C10_flyPAD.csv')
data_fp <- filter(data_fp, id == "sips_n")

data_fp$MS <- NA
data_fp$FS <- NA

data_fp$MS[grep('mated', data_fp$condition)] <- 'mated'
data_fp$MS[grep('virgin', data_fp$condition)] <- 'virgin'
data_fp$FS[grep('ff', data_fp$condition)] <- 'fed'
data_fp$FS[grep('10d', data_fp$condition)] <- 'deprived'

data_fp$FS <- factor(data_fp$FS, levels = c("fed", "deprived"))
data_fp$MS <- factor(data_fp$MS, levels = c("virgin", "mated"))
data_fp$spotFS <- paste0(data_fp$spot, ":", data_fp$FS)
data_fp$spotFS <-
  factor(
    data_fp$spotFS,
    levels = c(
      'sucrose:fed',
      'yeast:fed',
      'sucrose:deprived',
      'yeast:deprived'
    )
  )
data_fp$FSMS <- paste0(data_fp$FS, ":", data_fp$MS)
data_fp$FSMS <-
  factor(
    data_fp$FSMS,
    levels = c(
      'fed:virgin',
      'deprived:virgin',
      'fed:mated',
      'deprived:mated'
    )
  )

data_fp$spot <- factor(data_fp$spot, level = c('yeast', 'sucrose'))

stat_wilcox <- data_fp %>%
  group_by(spot) %>%
  wilcox_test(value ~ FSMS) %>%
  adjust_pvalue(method = 'holm') %>%
  add_significance()

plot_stats <- filter(stat_wilcox)
plot_stats <- plot_stats[-c(3, 4, 9, 10), ]
plot_stats$y.position <- c(14, 16.5, 19, 14,      1.7, 2, 2.3, 1.7)

labels_r <- c(
  'yeast:fed' = 'fully fed',
  'yeast:deprived' = 'protein deprived',
  'sucrose:fed' = 'fully fed',
  'sucrose:deprived' = 'protein deprived',
  
  'yeast' = 'protein',
  'sucrose' = 'sugar'
)


ggplot(filter(data_fp), aes(x = FSMS, y = value / 1000)) +
  geom_boxplot(
    aes(color = spot, fill = spot),
    width = .5,
    outlier.alpha = 1,
    outlier.size = .2,
    size = .3
  ) +
  stat_summary(
    color = "grey10",
    geom = "crossbar",
    width = .5,
    size = .2,
    fun.data = function(x) {
      return(c(
        y = median(x),
        ymin = median(x),
        ymax = median(x)
      ))
    }
  ) +
  scale_color_manual(values = c('#e4ba86ff', '#93d2ebff'), NULL) +
  scale_fill_manual(values = c('#e4ba86ff', '#93d2ebff'), NULL) +
  facet_grid(spot ~ ., scales = "free", labeller = as_labeller(labels_r)) +
  labs(y = 'sips on food (x1000)',
       fill = 'internal-state',
       color = 'internal-state',
       x = '') +
  theme_dahaniel2(base_size = 5, base_family = 'Arial') +
  theme(
    strip.text = element_text(size = 6),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(colour = "grey80"),
    axis.line.y = element_line(colour = "grey80"),
    legend.position = 'none'
  ) +
  coord_cartesian(clip = 'off') +
  stat_pvalue_manual(
    plot_stats,
    label = "{formatC(p, 2)}",
    size = stat_size,
    bracket.size = stat_bracket.size,
    tip.length = stat_tip.length,
    y.position = plot_stats$y.position
  )
ggsave(
  paste0(out_path, "/behavior-matrix.pdf"),
  width = 35,
  height = 40,
  dpi = plot_dpi,
  scale = plot_scale,
  plot = last_plot(),
  units = plot_units
)
write.csv(
  data_fp %>% select(spot, MS, FS, value) %>% arrange(spot),
  paste0(out_path, "/source-data_Fig1a.csv")
)
