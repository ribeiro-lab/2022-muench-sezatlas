# regions with sign. effects (wilcox) --------------------------------

data_peak <-
  peak_response(
    data,
    method = "mean",
    frames = 1:5,
    bgf = -5:0,
    sd_mp = SD_cutoff
  )
data_peak <- filter(data_peak, stimulus %ni% c('h2o', 'nostim'))

data_peak$MS <- factor(data_peak$MS, levels = c("virgin", "mated"))
data_peak$FS <- factor(data_peak$FS, levels = c("fed", "deprived"))
data_peak$MSFS <- paste0(data_peak$MS, ":", data_peak$FS)
data_peak$MSFS <-
  factor(
    data_peak$MSFS,
    levels = c(
      'virgin:fed',
      'virgin:deprived',
      'mated:fed',
      'mated:deprived'
    )
  )
data_peak$stimulus <-
  factor(data_peak$stimulus,
         levels = c('yeast', 'sucrose', 'nostim', 'h2o'))
data_peak$region <- factor(data_peak$region)

plot_data <- data_peak

#order
order_median <-
  plot_data %>% filter(stimulus == 'yeast', FS == 'fed', MS == 'mated') %>%
  group_by(region) %>%
  summarize(median = median(peak)) %>%
  arrange(median) %>%
  mutate(order_median = row_number())

order_delta <-
  plot_data %>% filter(stimulus == 'yeast', MS == 'mated') %>%
  group_by(region, FS) %>%
  summarize(median = median(peak)) %>%
  tidyr::spread(FS, median) %>%
  mutate(delta = deprived - fed) %>%
  arrange(delta) %>%
  ungroup %>%
  mutate(order_delta = dplyr::row_number()) %>%
  select(region, order_delta, delta)

plot_data <- merge(plot_data, order_median)
plot_data <- merge(plot_data, order_delta)

stat_wilcox <- plot_data %>%
  filter(MS == "mated") %>%
  group_by(region, stimulus) %>%
  wilcox_test(peak ~ FS) %>%
  add_significance() %>%
  add_xy_position(x = "region", dodge = 0.85) %>%
  filter(p.signif != 'ns')


p <-
  ggplot(data = filter(plot_data, MS == 'mated'), aes(x = reorder(region, order_delta), y = peak)) +
  geom_point(
    aes(color = FS),
    position = position_dodge(width = .85),
    alpha = .4,
    shape = 16,
    size = 2
  ) + #), size = 3.5, shape = 16, alpha = .6, cex = 5) +
  geom_boxplot(
    aes(color = FS),
    position = position_dodge(width = .85),
    size = .6,
    outlier.shape = 21,
    outlier.alpha = 1,
    outlier.fill = 'white',
    outlier.colour = 'grey60',
    alpha = 0
  ) +
  scale_color_manual(values = fs_colors, 'metabolic state') +
  scale_fill_manual(values = fs_colors, 'metabolic state') +
  facet_grid(. ~ stimulus) +
  labs(x = 'atlas region', y = lab_df) +
  theme_dahaniel2(base_size = 28) +
  theme(
    axis.text.x  = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = .5
    ),
    strip.text = element_text(size = 24),
    panel.grid.major.x = element_line(size = .3, color = 'grey80'),
    panel.grid.major.y = element_line(size = .2, color = 'grey80'),
    panel.border = element_blank(),
    legend.position = 'bottom'
  ) +
  coord_flip(ylim = c(-.4, 1.5))

p <-
  p + stat_pvalue_manual(
    stat_wilcox,
    label = "p.signif",
    x = 'region',
    size = 9,
    y.position = -.3,
    hjust = 0
  )


ggsave(
  paste0(out_path, "all_boxplots.pdf"),
  width = 150,
  height = 280,
  dpi = 300,
  scale = 3,
  units = 'mm',
  plot = p
)

write.csv(
  filter(plot_data, MS == 'mated', stimulus == 'yeast') %>%
    select(region, FS, peak) %>%
    arrange(FS, region),
  paste0(out_path, "/source-data_FigS8a.csv")
)

write.csv(
  filter(plot_data, MS == 'mated', stimulus == 'sucrose') %>%
    select(region, FS, peak) %>%
    arrange(FS, region),
  paste0(out_path, "/source-data_FigS8b.csv")
)



# Table S1 pooled fdr -----------------------------------------------------
pooled_data <- data_peak

stat_wilcox_y <- pooled_data %>%
  filter(MS == "mated", stimulus == 'yeast') %>%
  group_by(region, stimulus) %>%
  wilcox_test(peak ~ FS) %>%
  add_significance()

stat_wilcox_s <- pooled_data %>%
  filter(MS == "mated", stimulus == 'sucrose') %>%
  group_by(region, stimulus) %>%
  wilcox_test(peak ~ FS) %>%
  add_significance()


stat_wilcox_fdr_y <- pooled_data %>%
  filter(stimulus == 'yeast') %>%
  group_by(region, stimulus) %>%
  wilcox_test(peak ~ FS) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()


stat_wilcox_fdr_s <- pooled_data %>%
  filter(stimulus == 'sucrose') %>%
  group_by(region, stimulus) %>%
  wilcox_test(peak ~ FS) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()

write.csv(stat_wilcox_y,
          paste0(out_path, "stat_nocor_mated_yeast.csv"))
write.csv(stat_wilcox_s,
          paste0(out_path, "stat_nocor_mated_sucrose.csv"))
write.csv(stat_wilcox_fdr_y,
          paste0(out_path, "stat_fdr_pooled_yeast.csv"))
write.csv(stat_wilcox_fdr_s,
          paste0(out_path, "stat_fdr_pooled_sucrose.csv"))

delta_mated_y <-
  plot_data %>% filter(stimulus == 'yeast', MS == 'mated') %>%
  group_by(region, FS) %>%
  summarize(median = median(peak)) %>%
  tidyr::spread(FS, median) %>%
  mutate(delta = deprived - fed) %>%
  ungroup

delta_mated_s <-
  plot_data %>% filter(stimulus == 'sucrose', MS == 'mated') %>%
  group_by(region, FS) %>%
  summarize(median = median(peak)) %>%
  tidyr::spread(FS, median) %>%
  mutate(delta = deprived - fed) %>%
  ungroup

delta_pooled_y <- plot_data %>% filter(stimulus == 'yeast') %>%
  group_by(region, FS) %>%
  summarize(median = median(peak)) %>%
  tidyr::spread(FS, median) %>%
  mutate(delta = deprived - fed) %>%
  ungroup

delta_pooled_s <- plot_data %>% filter(stimulus == 'sucrose') %>%
  group_by(region, FS) %>%
  summarize(median = median(peak)) %>%
  tidyr::spread(FS, median) %>%
  mutate(delta = deprived - fed) %>%
  ungroup

write.csv(delta_mated_y, paste0(out_path, "delta_mated_yeast.csv"))
write.csv(delta_mated_s, paste0(out_path, "delta_mated_sucrose.csv"))
write.csv(delta_pooled_y, paste0(out_path, "delta_pooled_yeast.csv"))
write.csv(delta_pooled_s, paste0(out_path, "delta_pooled_sucrose.csv"))
